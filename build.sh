#!/bin/bash

# Delete existing pquant-conda
conda deactivate pquant-conda
conda remove --name pquant-conda --all

# build an R environment with base R
conda env create --quiet -f environment.yml && conda clean -a

# activate the environment
conda activate pquant-conda

# install R packages
conda install -c conda-forge r-devtools
R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("limma")

install.packages("ggplot2")

devtools::install_github("bartongroup/proteusLabelFree")
devtools::install_github("bartongroup/proteusTMT")
devtools::install_github("bartongroup/proteusSILAC")

devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
sessionInfo()
