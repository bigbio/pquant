## Data source
[RPXD015270.1-cell-lines](https://ftp.pride.ebi.ac.uk/pride/data/proteomes/proteogenomics/differential-expression/RPXD015270.1-cell-lines/expdesign/)
## MSstats
MSstats is an R-based/Bioconductor package for statistical relative quantification of peptides and proteins in mass spectrometry-based proteomic experiments. It is applicable to 
multiple types of sample preparation, including label-free workflows, workflows that use stable isotope labeled reference proteins and peptides, and work-flows that use fractionation. 
It is applicable to targeted Selected Reactin Monitoring(SRM), Data-Dependent Acquisiton(DDA or shotgun), and Data-Independent Acquisition(DIA or SWATH-MS).<br>
<br>
- After running either `RStudio_prePquant.R` or `prePquant.R`, we have all the ready data required by MSStats.
- Next to run `Rstudio_pquant.R`, we can see the data related **volcano plots, heatmap** and **QC plots** in `MSstats method` option.<br>
![](https://github.com/Douerww/pquantR/blob/main/img/MSstatsPage.png)

## Proteus
*Proteus* is an R package for downstream analysis of *MaxQuant* output. The input for *Proteus* is the evidence file. Evidence data are aggregated into peptides and then into 
proteins. *Proteus* offers many visualisation and data analysis tools both at peptide and protein level. In particular it allows simple differential expression using *limma*.<br>
<br>

- After running either `RStudio_prePquant.R` or `prePquant.R`, we could get a part of data required by Proteus.
- Next to run `Rstudio_pquant.R`, we can see the preprocessed data in `Proteus method` option.
- Follow the ShinyApp sidebar and annotate in the `Data` box，we can see the data related **normalization, mean variance relationship** and **protein clustering** in `Plots` box.<br>
![](https://github.com/Douerww/pquantR/blob/main/img/proteusPlotsPage.png)
