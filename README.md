# GEO_Data_Downloader
A highly configurable NCBI GEO database data download script that supports custom search keywords, download quantity, and number of concurrent threads. It allows personalized data retrieval and downloading solely through the configuration file without modifying the main script.

Before you run it, you need to set these configurations in download_config.py

///  
NCBI_EMAIL --> The email you registered with NCBI

NCBI_API_KEY --> Your NCBI API KEY

NCBI_KEYWORD --> GEO Search Keywords

MAX_GSE  --> Maximum number of search results

THREADS  --> Number of concurrent threads

OUT_ROOT --> Output root directory

///  

Run the script

Front-end operation：python Get_data_new.py

Run cluster in the background：nohup python Get_data.py > geo_download.log 2>&1 &

///  

Running result  
├── Downloaded_Data/       # Output root directory (can be customized in down_config.py as OUT_ROOT)  
│   ├── scRNA/                 # scRNA-seq data (automatically identified)    
│   │   ├── GSExxxxxx/   
│   │   │   ├── Metadata/      # Sample Metadata（full_metadata.csv）  
│   │   │   ├── Processed_Data/# Expression matrix（expression_matrix.txt）    
│   │   │   └── Supplementary_Raw/ # RAW.tar Original File    
│   ├── bulk/                  # bulk RNA-seq data  
│   ├── ATAC/                  # ATAC-seq data   
│   ├── ChIP/                  # ChIP-seq data  
│   ├── Unknown/               # Data of unrecognized technology type    
│   └── manifest.csv           # Download Status Record (Core of Resumable Download)    
├── config.py                  # Your custom configuration file  
├── Get_data_new.py            # Main Download Script    
└── geo_download.log           # Operation log (generated during background operation)    

///    

Attention!!!   

Since the GEO database can easily interrupt the download of large files, scRNA-seq downloads may sometimes fail. Please check the log file after running the script.   

Errors in determining the type of technology may occur due to misleading series titles or descriptions, so manual data verification is required after using this script.  

The script only supports retrieving GSE for organisms that are Homo sapiens. If you need to change the organism, please modify the function search_gse_list().    

///  
