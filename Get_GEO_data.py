#!/usr/bin/env python3
"""
GEO Self-Download Database Construction Script
- scRNA: Preferentially download only GSE*_RAW.tar (reduces fragmented downloads for each GSM)
- Supports resume from breakpoints (Range), timeout + retry, and atomic renaming of temporary files
- manifest.csv is used to record download status, making it easy to resume after interruption
- Thread concurrency is controllable
"""

import os
import time
import gzip
import re
import socket
import csv
from pathlib import Path
from random import uniform
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
import pandas as pd
import numpy as np
from Bio import Entrez
from tqdm import tqdm

# ===================== User Settings =====================
try:
    from config import NCBI_EMAIL, NCBI_API_KEY, NCBI_KEYWORD, MAX_GSE, THREADS, OUT_ROOT
    Entrez.email = NCBI_EMAIL
    Entrez.api_key = NCBI_API_KEY
    KEYWORD = NCBI_KEYWORD
    MAX_GSE = MAX_GSE
    THREADS = THREADS
    OUT_ROOT = OUT_ROOT
except ImportError:
    Entrez.email = "your_email@example.com"
    Entrez.api_key = ""
    KEYWORD = ""
    MAX_GSE = 1
    THREADS = 1
    OUT_ROOT = Path("GEO_Database_Raw")
    print("Please create a config.py file and fill in NCBI_EMAIL and NCBI_API_KEY!")


# Technical Keyword Definition
TECH_KEYWORDS = {
    "scRNA": ["single cell", "single-cell", "scrna", "10x", "chromium", "drop-seq", "cite-seq"],
    "ATAC": ["atac", "accessibility", "chromatin"],
    "ChIP": ["chip-seq", "chipseq", "histone", "tf binding"],
    "bulk": ["rna-seq", "transcriptome", "expression profiling"]
}

# Expected file extension classification
FILE_TYPES = {
    "counts": ["count", "counts", "raw", "read_count"],
    "peaks": ["peak", "narrowpeak", "broadpeak", "bed"],
    "fragments": ["fragment", "fragments"],
    "annotation": ["annotation", "gtf", "gff", "gene_info"]
}

# Internet settings
socket.setdefaulttimeout(120)
HEADERS = {"User-Agent": "Mozilla/5.0 (compatible; GEO-Downloader/1.0)"}

# Initialize Directory & manifest
for tech in list(TECH_KEYWORDS.keys()) + ["Unknown"]:
    (OUT_ROOT / tech).mkdir(parents=True, exist_ok=True)

MANIFEST = OUT_ROOT / "manifest.csv"
MANIFEST_FIELDS = [
    "GSE", "Technology", "SeriesMatrix", "SeriesMatrix_Path",
    "RAW_downloaded", "RAW_path", "Metadata_path", "Processed_files", "Notes"
]

# Read existing manifest (for resumable transfers)
def load_manifest():
    if MANIFEST.exists():
        try:
            return pd.read_csv(MANIFEST)
        except:
            return pd.DataFrame(columns=MANIFEST_FIELDS)
    else:
        return pd.DataFrame(columns=MANIFEST_FIELDS)


def save_manifest_row(row):
    """Add or update a record about a specific GSE in the manifest"""
    df = load_manifest()
    gse = row.get("GSE")

    clean_row = {}
    for k, v in row.items():
        clean_row[k] = np.nan if v == '' else v

    if gse in df.get("GSE", []).values:
        df.loc[df["GSE"] == gse, list(clean_row.keys())] = list(clean_row.values())
    else:
        df = pd.concat([df, pd.DataFrame([clean_row])], ignore_index=True)

    df.to_csv(MANIFEST, index=False)

# ------------------- Function -------------------

def get_gse_prefix(gse):
    """GSE123456 -> GSE123nnn"""
    gse_num = gse[3:]
    if len(gse_num) <= 3:
        return "GSEnnn"
    return f"GSE{gse_num[:-3]}nnn"

def infer_tech(info):
    text = (info.get("title", "") + " " + info.get("Type", "") + " " + info.get("summary", "")).lower()
    for tech, kws in TECH_KEYWORDS.items():
        if any(kw in text for kw in kws):
            return tech
    return "Unknown"

# Robust download with resume (Range), retries, atomic rename
def robust_download(url, dest_path, max_retries=6, chunk_size=4*1024*1024, min_valid_size=100):
    dest = Path(dest_path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")

    # NCBI ftp -> https
    if str(url).startswith("ftp://"):
        url = url.replace("ftp://", "https://")

    session = requests.Session()
    headers = HEADERS.copy()
    attempt = 0
    backoff = 5

    while attempt < max_retries:
        attempt += 1
        try:
            resume_byte_pos = tmp.stat().st_size if tmp.exists() else 0
            if resume_byte_pos > 0:
                headers["Range"] = f"bytes={resume_byte_pos}-"
            else:
                headers.pop("Range", None)

            with session.get(url, stream=True, headers=headers, timeout=(30, 300)) as r:
                if r.status_code in (403, 404):
                    # Not found or forbidden: no reason to retry a lot
                    # but may be transient - allow a couple retries
                    if attempt >= max_retries:
                        return False
                # Accept 200 or 206 for resume
                if r.status_code not in (200, 206):
                    # may be redirect to HTML error; treat as fail & retry
                    raise RuntimeError(f"HTTP {r.status_code}")

                total = None
                content_range = r.headers.get("Content-Range")
                if content_range:
                    # Content-Range: bytes 100-999/1000
                    try:
                        total = int(content_range.split("/")[-1])
                    except:
                        total = None
                else:
                    cl = r.headers.get("Content-Length")
                    if cl:
                        try:
                            total = int(cl) + resume_byte_pos
                        except:
                            total = None

                mode = "ab" if resume_byte_pos > 0 else "wb"
                with open(tmp, mode) as f:
                    # stream in large chunks
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                # finished this attempt's streaming
            final_size = tmp.stat().st_size if tmp.exists() else 0
            if final_size == 0 or (min_valid_size and final_size < min_valid_size):
                # small or zero file -> treat as fail
                raise RuntimeError(f"Downloaded size too small: {final_size}")
            # success -> rename
            tmp.replace(dest)
            return True
        except Exception as e:
            # on error, wait and retry
            # print minimal message to avoid flooding
            print(f"[Retry {attempt}/{max_retries}] {Path(dest).name}: {e}")
            time.sleep(backoff)
            backoff = min(backoff * 2, 300)
            continue
    # all attempts failed
    if tmp.exists():
        try:
            tmp.unlink()
        except:
            pass
    return False

# ------------------- series_matrix Analysis -------------------

def parse_series_matrix(gse_id, matrix_path, output_dir):
    output_dir = Path(output_dir)
    meta_out = output_dir / "Metadata"
    data_out = output_dir / "Processed_Data"
    supp_out = output_dir / "Supplementary_Raw"
    for d in [meta_out, data_out, supp_out]:
        d.mkdir(parents=True, exist_ok=True)

    ftp_urls = []
    gsm_ids = []
    meta_dict = {}
    matrix_start = -1
    matrix_end = -1

    try:
        with gzip.open(matrix_path, "rt", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
    except Exception as e:
        print(f"[Parse Error] {gse_id} cannot read series matrix: {e}")
        return [], None, None

    for idx, raw in enumerate(lines):
        line = raw.strip()
        # collect urls
        for u in re.findall(r'(https?://\S+|ftp://\S+)', line):
            u = u.strip().strip('"').rstrip(',')
            if "series_matrix.txt.gz" not in u:
                ftp_urls.append(u)
        # sample ids
        if line.startswith("!Sample_geo_accession"):
            parts = line.split("\t")
            gsm_ids = [p.strip().strip('"') for p in parts[1:]]
        elif line.startswith("!Sample_"):
            parts = line.split("\t")
            raw_key = parts[0].strip().strip('"').replace("!Sample_", "")
            vals = [p.strip().strip('"') for p in parts[1:]]
            if gsm_ids and len(vals) != len(gsm_ids):
                vals += [""] * (len(gsm_ids) - len(vals))
                vals = vals[:len(gsm_ids)]
            key = raw_key
            cnt = 1
            while key in meta_dict:
                cnt += 1
                key = f"{raw_key}_{cnt}"
            meta_dict[key] = vals
        elif line.startswith("!series_matrix_table_begin"):
            matrix_start = idx + 1
        elif line.startswith("!series_matrix_table_end"):
            matrix_end = idx
            break

    # write metadata CSV
    metadata_path = None
    if gsm_ids:
        df_meta = pd.DataFrame({"GSM": gsm_ids})
        for k, v in meta_dict.items():
            if len(v) == len(gsm_ids):
                df_meta[k] = v
        metadata_path = meta_out / "full_metadata.csv"
        df_meta.to_csv(metadata_path, index=False)

    # write expression matrix if present
    expression_path = None
    if matrix_start > 0 and matrix_end > matrix_start:
        header = lines[matrix_start].split("\t")
        if len(header) > 1:
            tmp_mat = Path(matrix_path).with_suffix(".tmp_matrix.txt")
            with open(tmp_mat, "w", encoding="utf-8") as tf:
                tf.writelines(lines[matrix_start:matrix_end])
            try:
                df_expr = pd.read_csv(tmp_mat, sep="\t", index_col=0, low_memory=False)
                df_expr.columns = [c.strip().strip('"') for c in df_expr.columns]
                expression_path = data_out / "expression_matrix.txt"
                df_expr.to_csv(expression_path, sep="\t")
            except Exception as e:
                print(f"[Warn] {gse_id} expression parse failed: {e}")
            finally:
                try:
                    tmp_mat.unlink()
                except:
                    pass

    ftp_urls = list(dict.fromkeys(ftp_urls))  # dedupe while preserving order
    return ftp_urls, metadata_path, expression_path

# ------------------- 处理单个 GSE -------------------

def process_gse(info):
    gse = info.get("Accession")
    print(f"\n[PROCESS] {gse} ...")
    tech = infer_tech(info)
    base_dir = OUT_ROOT / tech / gse
    meta_dir = base_dir / "Metadata"
    data_dir = base_dir / "Processed_Data"
    supp_dir = base_dir / "Supplementary_Raw"
    for d in [meta_dir, data_dir, supp_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # record initial manifest row
    manifest_row = {
        "GSE": gse,
        "Technology": tech,
        "SeriesMatrix": False,
        "SeriesMatrix_Path": "",
        "RAW_downloaded": False,
        "RAW_path": "",
        "Metadata_path": "",
        "Processed_files": "",
        "Notes": ""
    }
    save_manifest_row(manifest_row)

    prefix = get_gse_prefix(gse)
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/{gse}_series_matrix.txt.gz"
    matrix_local = base_dir / f"{gse}_series_matrix.txt.gz"
    if not matrix_local.exists():
        ok = robust_download(matrix_url, matrix_local)
    else:
        ok = True

    if ok and matrix_local.exists():
        manifest_row["SeriesMatrix"] = True
        manifest_row["SeriesMatrix_Path"] = str(matrix_local)
        save_manifest_row(manifest_row)
    else:
        manifest_row["Notes"] += "series_matrix_missing; "
        save_manifest_row(manifest_row)

    ftp_urls, metadata_path, expression_path = [], None, None
    if matrix_local.exists():
        ftp_urls, metadata_path, expression_path = parse_series_matrix(gse, matrix_local, base_dir)
        if metadata_path:
            manifest_row["Metadata_path"] = str(metadata_path)
        if expression_path:
            manifest_row["Processed_files"] = "expression_matrix"
        save_manifest_row(manifest_row)

    try:
        sra_meta = {}
        # lightweight: use Entrez.esearch + esummary on sra (may be limited)
        term = f"{gse}[Accession]"
        h = Entrez.esearch(db="sra", term=term, retmax=200)
        ids = Entrez.read(h).get("IdList", [])
        h.close()
        if ids:
            h2 = Entrez.esummary(db="sra", id=",".join(ids))
            recs = Entrez.read(h2)
            h2.close()
            # Try extract runs and layout
            for rec in recs:
                try:
                    expxml = rec.get("ExpXml", "")
                    gsmm = re.search(r"(GSM\d+)", expxml)
                    if gsmm:
                        gsm = gsmm.group(1)
                        layout = "Unknown"
                        if 'LibraryLayout="PAIRED"' in expxml: layout = "PAIRED"
                        if 'LibraryLayout="SINGLE"' in expxml: layout = "SINGLE"
                        platform = rec.get("Platform", "Unknown")
                        sra_meta[gsm] = {"layout": layout, "platform": platform, "raw": rec.get("Runs", "")}
                except:
                    continue
    except Exception:
        sra_meta = {}

    if tech == "scRNA":
        # download RAW.tar only
        raw_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/{gse}_RAW.tar"
        raw_local = supp_dir / f"{gse}_RAW.tar"
        if not raw_local.exists():
            ok_raw = robust_download(raw_url, raw_local, chunk_size=1*1024*1024, max_retries=10)
        else:
            ok_raw = True
        if ok_raw and raw_local.exists():
            manifest_row["RAW_downloaded"] = True
            manifest_row["RAW_path"] = str(raw_local)
            manifest_row["Notes"] += "scRNA_RAW_saved; "
        else:
            manifest_row["Notes"] += "scRNA_RAW_missing; "
        save_manifest_row(manifest_row)
        # Do not iterate per-GSM downloads for scRNA to avoid fragmentation
        return f"{gse} [{tech}] done"

    processed_files = []
    if ftp_urls:
        supp_dir.mkdir(parents=True, exist_ok=True)
        for url in ftp_urls:
            fname = url.split("/")[-1]
            # filter routine files
            low = fname.lower()
            if any(x in low for x in ['.bmp', '.jpg', 'readme', 'filelist']):
                continue
            # try download important file types
            if any(x in low for x in ['count','matrix','mtx','tsv','csv','xls','bed','narrowPeak','broadPeak','fragments','fragment']):
                dest = supp_dir / fname
                if robust_download(url, dest, max_retries=4):
                    processed_files.append(str(dest))
    if tech in ["ATAC", "ChIP", "bulk"]:
        raw_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/{gse}_RAW.tar"
        raw_local = supp_dir / f"{gse}_RAW.tar"
        if not raw_local.exists():
            ok_raw = robust_download(raw_url, raw_local, max_retries=4)
        else:
            ok_raw = True
        if ok_raw and raw_local.exists():
            manifest_row["RAW_downloaded"] = True
            manifest_row["RAW_path"] = str(raw_local)
            processed_files.append(str(raw_local))

    if processed_files:
        manifest_row["Processed_files"] = ";".join(processed_files)
    save_manifest_row(manifest_row)
    return f"{gse} [{tech}] done"

# ------------------- Main Process -------------------
def search_gse_list(keyword, max_gse=1000):
    print(f"[SEARCH] searching: {keyword}")
    query = f"({keyword}) AND \"Homo sapiens\"[Organism] AND GSE[Entry Type]"

    try:
        # 第一次搜索仅为了获取总数
        h = Entrez.esearch(db="gds", term=query, retmax=0)
        total = int(Entrez.read(h)["Count"])
        h.close()

        limit = min(total, max_gse)
        print(f"[SEARCH] found {total}, processing first {limit}")
    except Exception as e:
        print("[ERROR] Entrez search failed:", e)
        return []

    results = []
    retstart = 0
    batch = 100

    while retstart < limit:
        try:
            # 1. 获取 ID 列表 (不需要 usehistory)
            h = Entrez.esearch(db="gds", term=query, retstart=retstart, retmax=min(batch, limit - retstart))
            res = Entrez.read(h)
            h.close()

            ids = res["IdList"]
            if not ids:
                break

            # 2. 修复点：直接使用 ID 列表请求 summary，而不是依赖 WebEnv 的偏移量
            # 这样确保获取的就是本轮循环查到的那些 ID
            h2 = Entrez.esummary(db="gds", id=",".join(ids))
            recs = Entrez.read(h2)
            h2.close()

            for r in recs:
                results.append({
                    "Accession": r.get("Accession"),
                    "Type": r.get("gdsType", ""),
                    "n_samples": int(r.get("n_samples", 0)),
                    "title": r.get("title", ""),
                    "summary": r.get("summary", "")
                })

            retstart += len(ids)
            print(f"  > Fetched batch {retstart}/{limit}")  # 增加进度提示
            time.sleep(0.4)

        except Exception as e:
            print(f"[WARN] batch fetch failed at retstart={retstart}: {e}")
            time.sleep(5)
            # 即使失败也尝试跳过当前 batch，避免死循环
            retstart += batch
            continue

    return results

def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print("Starting GEO downloader (improved).")
    existing_manifest = load_manifest()
    gse_list = search_gse_list(KEYWORD, max_gse=MAX_GSE)
    if not gse_list:
        print("No GSEs found, exit.")
        return

    final_list = []
    processed = set()
    for item in gse_list:
        gse = item["Accession"]
        # get subseries
        try:
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief"
            r = requests.get(url, headers=HEADERS, timeout=10)
            if r.status_code == 200 and "SuperSeries of:" in r.text:
                subs = []
                for line in r.text.splitlines():
                    if "SuperSeries of:" in line:
                        subs.extend(re.findall(r"GSE\d+", line))
                for s in sorted(set(subs)):
                    if s not in processed:
                        final_list.append({"Accession": s, "Type":"Unknown", "n_samples":0, "title": item.get("title",""), "summary": ""})
                        processed.add(s)
                processed.add(gse)
                continue
        except:
            pass
        if gse not in processed:
            final_list.append(item)
            processed.add(gse)

    print(f"[INFO] total to process: {len(final_list)}")

    # 3. parallel processing (manifest allows resume)
    with ThreadPoolExecutor(max_workers=THREADS) as exe:
        futures = {exe.submit(process_gse, info): info for info in final_list}
        for f in tqdm(as_completed(futures), total=len(futures), desc="Downloading"):
            info = futures[f]
            try:
                res = f.result()
                print("  >", res)
            except Exception as e:
                print(f"[ERROR] {info.get('Accession')} failed: {e}")

    print("All done. Manifest:", MANIFEST)

if __name__ == "__main__":
    main()
