Bootstrap: docker
From: rocker/tidyverse:4.3.1

%labels
    Version v0.0.1

%environment
    # Persistent environment variables
    export R_LIBS_USER=/usr/local/lib/R/site-library
    export R_LIBS=/usr/local/lib/R/site-library
    export PATH=/usr/local/bin:$PATH
    export R_ENVIRON_USER=/usr/local/lib/R/etc/Renviron.site

%post
    # BiocManager to install packages
    R -e "library("BiocManager"); BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'data.table', 'MungeSumstats', 'rtracklayer', 'SeqArray', 'SNPRelate'), dependencies=TRUE)"
    R -e "install.packages(c('argparse', 'coloc', 'FactoMineR', 'foreach', 'future.apply', 'ggforce', 'ggrepel', 'httr', 'parallel', 'patchwork', 'purr', 'stringi', 'susieR'), dependencies=TRUE)"
    R -e "devtools::install_github('heatherjzhou/PCAForQTL')"

%runscript
    echo "Running TOPCHeF R environment..."
    exec R "$@"

%help
    This is a definition file for all R/4.3.1 TOPCHeF analyses.