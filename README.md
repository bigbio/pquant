# pquantR
pquantR is a R shiny application which is designed to do the downstream analysis of proteomics dataset, currently these figures are included: Heatmap, Volcano Plot, QC plot.<br>
Because the application is in developing, the test figures are drawn by test R packages separately now.<br>
## Test procedure
1. Preparing the test data for proteus.<br>
* Download out.csv and out.mzTab from [this page](http://ftp.pride.ebi.ac.uk/pride/data/proteomes/RPXD012431.1/proteomics_lfq/).
* Run get_proteus_evidence.py, it will combine out.csv and out.mzTab to obtain the evidence.csv, as input of proteus. (It still have some bugs)
* The metadata.csv (for proteus) is manual annotated in “datasets”. We will code a python file to get it in the future.

2. Use msstat to draw the figures.
* Run test_MSstats.R.<br>
`tips: You should put the R scripts and the test data in a same path, if you want to put it in a different path, please modify the code that reads the path`<br>
test_MSstats.R can draw many plots such as Heatmap, Volcano Plot, QC plot, Condition plot, Comparison plot and so on.<br>
![](https://github.com/Douerww/pquantR/blob/main/image/MSstats_output.png)
3. Use proteus to draw the figures.
* Run test_proteus.R<br>
`tips: You should put the R scripts and the test data in a same path, if you want to put it in a different path, please modify the code that reads the path`<br>
test_proteus.R can create a peptide dataset, create protein dataset, aggregate protein and draw a volcano plot. The path where you store this R script.<br>
![](https://github.com/Douerww/pquantR/blob/main/image/proteus_output.png)
