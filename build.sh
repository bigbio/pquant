#!/bin/bash

# Delete existing pquant-conda
conda remove --name pquant-conda --all

# build an R environment with base R
conda env create --quiet -f environment.yml && conda clean -a

# activate the environment
conda activate pquant-conda

# Active R environment.
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
R -e "BiocManager::install()"
R -e "BiocManager::install('limma')"
R -e "install.packages('ggplot2')"
R -e "devtools::install_github('bartongroup/proteusLabelFree')"
R -e "devtools::install_github('bartongroup/proteusTMT')"
R -e "devtools::install_github('bartongroup/proteusSILAC')"
R -e "devtools::install_github('bartongroup/Proteus')"
R -r "sessionInfo()"
