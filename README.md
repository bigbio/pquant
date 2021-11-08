# pquant

pquantR is a python and R package to perform downstream analysis of proteomicsLFQ quantitative data. It also included R shiny application which is designed to do the downstream analysis of proteomics dataset, currently these figures are included: Heatmap, Volcano Plot, QC plot.<br>

Because the application is in developing, the test figures are drawn by test R packages separately now.<br>

## Installation of environment

The `build.sh` bash script use conda and mamba to install the environment to install all the dependencies and packahes in python and R.

Prerequisites:

- conda (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- mamba (https://github.com/mamba-org/mamba)

After the installation of conda and mamba you can run the following script:

```bash
$> source build.sh
```
## sample datasets

- UPS1:
  - Configuration run from proteomicsLFQ  (unique peptides, star align, stricter pep and protein FDR) - [unique-peptides-star-align-stricter-pep-protein-FDR](https://ftp.pride.ebi.ac.uk/pride/data/proteomes/ups1/unique-peptides-star-align-stricter-pep-protein-FDR/proteomics_lfq/)


## Shiny application

1. Preparing  for shiny app.<br>
* Download `pquantR` folder from [this page](https://github.com/Douerww/pquantR/tree/main/pquant) and put it in a suitable working path.
2. Run `app.R` in the folder, and you could run it in the following ways.<br>
* Upload the data, choose the parameters of MSstats then submit.<br>
3. We could see visualization of processed data and differentially abundant proteins.<br>
![](https://github.com/Douerww/pquantR/blob/main/img/home.png)
![](https://github.com/Douerww/pquantR/blob/main/img/data-1.png)
![](https://github.com/Douerww/pquantR/blob/main/img/data-2.png)
![](https://github.com/Douerww/pquantR/blob/main/img/data-3.png)
![](https://github.com/Douerww/pquantR/blob/main/img/msstats-1.png)
![](https://github.com/Douerww/pquantR/blob/main/img/msstats-2.png)
![](https://github.com/Douerww/pquantR/blob/main/img/dynamic-1.png)
![](https://github.com/Douerww/pquantR/blob/main/img/dynamic-2.png)


## Todo list
1. Add download function, including PDF,PNG format, etc.<br>
2. Other interactive plots.<br>
3. Protein accession -> gene accession.<br>
4. Portein accessions with expression values and connect to REACTOME.<br>
