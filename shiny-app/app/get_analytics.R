#! /usr/bin/R

getAnalytics <- function(data){

#print(getwd())
setwd('../data/')

#default option
DDA2009.proposed <- dataProcess(raw = data,
                                normalization = 'equalizeMedians',
                                summaryMethod = 'TMP',
                                censoredInt = "NA",
                                cutoffCensored = "minFeature",
                                MBimpute = TRUE,
                                maxQuantileforCensored=0.999)

# No imputation
DDA2009.TMP <- dataProcess(raw = data,
                           normalization = 'equalizeMedians',
                           summaryMethod = 'TMP',
                           censoredInt = NULL, MBimpute=FALSE)


# Automatically create the manually created matrix in MSStats, User Manual P23
len <- length(levels(DDA2009.TMP$ProcessedData$GROUP_ORIGINAL))

ourMatrix <- matrix(c(0:0),nrow=len,ncol=len)
diag(ourMatrix) = -1
for(i in 1:len-1){
  ourMatrix[i,i+1] = 1
}
ourMatrix[len,1] = 1

ourCondition <- unique(data$Condition)
len2 <- length(ourCondition)
tmp <- matrix(ourCondition, nr=len2, nc=1)
name <- matrix(nr=len2, nc=1)
for(i in 1:len2-1){
  name[i,1] <- sprintf('%s-%s', tmp[i+1,1], tmp[i,1])
}
name[len2,1] <- sprintf('%s-%s', tmp[1,1], tmp[len2,1])

row.names(ourMatrix) <- name
#----------end-----------

DDA2009.comparisons <- groupComparison(contrast.matrix = ourMatrix,
                                       data = DDA2009.proposed)

write.csv(DDA2009.comparisons$ComparisonResult, file="MSstats_output.csv")


####
#data = DDA2009.comparisons$ComparisonResult

name_data <- read.table('g_g_name.txt', sep = ',')


for(i in 1:(length(name_data)/2)){
  data <- read.csv("MSstats_output.csv")
  
  temp.name = name_data[,i]
  data <- data[which(data$Label %in% temp.name), ]
  
  if (any(is.na(data))) {
      matrix_final = data.frame((Protein = data$Protein))
      colnames(matrix_final)[1] <- 'Protein'
      break
  }
}



for(i in 1:(length(name_data)/2)){
  data <- read.csv("MSstats_output.csv")
  
  temp.name = name_data[,i]
  data <- data[which(data$Label %in% temp.name), ]
  
  if (any(is.na(data))) {
      num = 1
      
      name_p = paste(name_data[,i+4], '.p-value', sep = '')
      name_log2fc = paste(name_data[,i+4], '.log2foldchange', sep = '')
      
      matrix_tmp = data.frame(pp = data$pvalue, llog2 = data$log2FC)
      colnames(matrix_tmp)[1] <- name_p
      colnames(matrix_tmp)[2] <- name_log2fc
      
      matrix_final = cbind(matrix_final, matrix_tmp)
  }
  else{
      next
  }
}


### Delete lines containing 'NA'
matrix_final = na.omit(matrix_final)

write.csv(matrix_final, file = 'NAME_analytics.csv', quote = F, row.names = F)

}


