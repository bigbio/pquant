#!/bin/bash

# Delete existing pquant-conda
conda remove --name pquant-conda --all

# build an R environment with base R
mamba env create --file environment.yml && conda clean -a -y

# activate the environment
conda activate pquant-conda

## Install pip dependencies
pip install rpy2
pip install numpy

# Active R environment.
R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
R -e "BiocManager::install()"
R -e "devtools::install_github('bartongroup/proteusLabelFree')"
R -e "devtools::install_github('bartongroup/proteusTMT')"
R -e "devtools::install_github('bartongroup/proteusSILAC')"
R -e "devtools::install_github('bartongroup/Proteus')"
R -e "devtools::install_github('Vitek-Lab/MSstats') "
R -e "install.packages('reticulate')"
R -r "sessionInfo()"
