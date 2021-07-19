#! /usr/bin/R


getPlot <- function(fileData, flag, selector, DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons){

# use type="QCplot" with all proteins
# change the upper limit of y-axis=35
# set up the size of pdf
if (flag == 'qc'){
  dataProcessPlots(data = DDA2009.proposed, type="QCplot",which.Protein=selector,
                 ylimDown=0, ylimUp=35,width=5, height=5, address=FALSE)
}

# volcanoPlots
if (flag == 'volcano'){
  groupComparisonPlots(data = DDA2009.comparisons$ComparisonResult, type = 'VolcanoPlot',
                       width=5, height=5, address=FALSE, which.Comparison=selector)
}

# Heatmaps
if (flag == 'heat'){
  heatmap <- read.csv('./pheatmap_input.csv', row.names = 1)
  
  if (nrow(heatmap) > 50) {pheatmap(heatmap, show_rownames = F)}
  else {pheatmap(heatmap)}
}


}