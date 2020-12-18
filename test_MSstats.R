#! /usr/bin/R

library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
setwd('D:/dataset/R downstream analysis/MSstats/output')  # user your own path
data <- read.csv("D:/dataset/R downstream analysis/MSstats/ourdata/out.csv") # from RPXD012431.1

#default option
DDA2009.proposed <- dataProcess(raw = data,
                                normalization = 'equalizeMedians',
                                summaryMethod = 'TMP',
                                censoredInt = "NA",
                                cutoffCensored = "minFeature",
                                MBimpute = TRUE,
                                maxQuantileforCensored=0.999)

# use type="QCplot" with all proteins
# change the upper limit of y-axis=35
# set up the size of pdf
dataProcessPlots(data = DDA2009.proposed, type="QCplot", ylimUp=35,
                 width=5, height=5)


# use type="QCplot" for 1st protein only
# change the upper limit of y-axis=35
# set up the size of pdf
dataProcessPlots(data = DDA2009.proposed, type="QCplot", which.Protein=1,
                 ylimUp=35, width=5, height=5)


# use type="Profileplot"
dataProcessPlots(data = DDA2009.proposed, type="Profileplot", ylimUp=35,
                 featureName="NA", width=5, height=5, address="DDA2009_proposed_")

# use type="Conditionplot"
dataProcessPlots(data = DDA2009.proposed, type="Conditionplot",
                 width=5, height=5, address="DDA2009_proposed_")

# No imputation
DDA2009.TMP <- dataProcess(raw = DDARawData,
                           normalization = 'equalizeMedians',
                           summaryMethod = 'TMP',
                           censoredInt = NULL, MBimpute=FALSE)

levels(DDA2009.TMP$ProcessedData$GROUP_ORIGINAL)

# The test data has 4 conditions
# The following comparison is rewritten to suit your individual needs, and next is just an example
comparison1 <- matrix(c(-1,1,0,0),nrow=1)
comparison2 <- matrix(c(0,-1,1,0),nrow=1)
comparison3 <- matrix(c(0,0,-1,1),nrow=1)
comparison4 <- matrix(c(1,0,0,-1),nrow=1)
comparison<-rbind(comparison1,comparison2,comparison3,comparison4)
row.names(comparison) <- c("C2-C1","C3-C2","C4-C3","C1-C4")
DDA2009.comparisons <- groupComparison(contrast.matrix = comparison, data = DDA2009.proposed)


# normal quantile-quantile plots
modelBasedQCPlots(data=DDA2009.comparisons, type="QQPlots",
                  width=5, height=5, address="DDA2009_proposed_")
# residual plots
modelBasedQCPlots(data=DDA2009.comparisons, type="ResidualPlots",
                  width=5, height=5, address="DDA2009_proposed_")

# volcanoPlots
groupComparisonPlots(data = DDA2009.comparisons$ComparisonResult, type = 'VolcanoPlot',
                     width=5, height=5, address="DDA2009_proposed_")

# Heatmaps
groupComparisonPlots(data = DDA2009.comparisons$ComparisonResult, type = 'Heatmap')

# ComparisonPlot
groupComparisonPlots(data=DDA2009.comparisons$ComparisonResult, type="ComparisonPlot",
                     width=5, height=5, address="DDA2009_proposed_")