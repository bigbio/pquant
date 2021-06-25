library('MSstats', warn.conflicts = F, quietly = T, verbose = F)
library(reticulate)

### get a "/pquant/data" path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./data/")
fileData <- read.csv('out_msstats.csv')

DDA2009.proposed <- dataProcess(raw = fileData,
                                normalization = 'equalizeMedians',
                                summaryMethod = 'TMP',
                                censoredInt = "NA",
                                cutoffCensored = "minFeature",
                                MBimpute = TRUE,
                                maxQuantileforCensored=0.999)


DDA2009.TMP <- dataProcess(raw = fileData,
                           normalization = 'equalizeMedians',
                           summaryMethod = 'TMP',
                           censoredInt = NULL, MBimpute=FALSE)


# Automatically create the manually created matrix in MSstats, user manual p23
len <- length(levels(DDA2009.TMP$ProcessedData$GROUP_ORIGINAL))

ourMatrix <- matrix(c(0:0),nrow=len,ncol=len)
diag(ourMatrix) = -1
for(i in 1:len-1){
  ourMatrix[i,i+1] = 1
}
ourMatrix[len,1] = 1

ourCondition <- unique(fileData$Condition)
len2 <- length(ourCondition)
tmp <- matrix(ourCondition, nr=len2, nc=1)
name <- matrix(nr=len2, nc=1)
for(i in 1:len2-1){
  name[i,1] <- sprintf('%s-%s', tmp[i+1,1], tmp[i,1])
}
name[len2,1] <- sprintf('%s-%s', tmp[1,1], tmp[len2,1])

row.names(ourMatrix) <- name
#----------End of creation-----------

DDA2009.comparisons <- groupComparison(contrast.matrix = ourMatrix,
                                       data = DDA2009.proposed)

save(DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons, file = "preShiny.RData")

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output.csv")

#! /usr/bin/python
#conda_install(packages = 'pandas') # If you are using it for the first time, you need to install the pandas package

py_run_file('../py/MSstatas to pheatmap.py')

py_run_file('../py/get_proteus_evidence.py')

proteus_evidence_file <- read.csv('out_proteus.csv', row.names = NULL)
# Prevent the occurrence of ''(null value)
proteus_evidence_file <- unique(proteus_evidence_file)
proteus_evidence_file[proteus_evidence_file == ''] <- NA
proteus_evidence_file <- na.omit(proteus_evidence_file)
write.csv(proteus_evidence_file, 'out_proteus.csv', row.names = F)

py_run_file('../py/generate_configuration_xml.py')
