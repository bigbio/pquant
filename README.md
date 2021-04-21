# pquantR
pquantR is a R shiny application which is designed to do the downstream analysis of proteomics dataset, currently these figures are included: Heatmap, Volcano Plot, QC plot.<br>
Because the application is in developing, the test figures are drawn by test R packages separately now.<br>
## Test procedure
1. Preparing  for shiny app.<br>
* Download `ups1-uniprot-comet-proteomics_lfq.zip` from [this page](https://github.com/bigbio/pquant/issues/7), we need `out_msstats.csv` in this zip.<br>
* Download `shiny-app` folder from [this page](https://github.com/Douerww/pquantR/tree/main/shiny-app) and put it in a suitable working path.
2. Run `app.R` in the shiny-app folder, and you could run it in the following ways.<br>
* Use RStudio to open app.R file, and click the `Run App` button.<br>
* Or enter the following command in RStudio to run the shiny app:<br>
```r
setwd(getwd())  # The path where you store this R script
shiny::runApp('./shiny-app')
```
3. Click the `Browse` button to upload ‘out. CSV ‘file, then you can view the contents of the file, the resulting volcano map, heat map, and QC map through four different buttons in the shiny app. And you could use MSstats method or proteus method.<br>
![](https://github.com/Douerww/pquantR/blob/main/img/homePage.png)
![](https://github.com/Douerww/pquantR/blob/main/img/defaultMethodPage.png)
![](https://github.com/Douerww/pquantR/blob/main/img/proteusDataPage.png)
![](https://github.com/Douerww/pquantR/blob/main/img/proteusPlotsPage.png)




