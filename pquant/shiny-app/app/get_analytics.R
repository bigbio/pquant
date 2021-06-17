#! /usr/bin/R

getAnalytics <- function(data, DDA2009.proposed, DDA2009.TMP, DDA2009.comparisons){



####
#data = DDA2009.comparisons$ComparisonResult

name_data <- read.table('data/g_g_name.txt', sep = ',')


for(i in 1:(length(name_data)/2)){
  data <- read.csv("data/MSstats_output.csv")
  
  temp.name = name_data[,i]
  data <- data[which(data$Label %in% temp.name), ]
  
  if (any(is.na(data))) {
      matrix_final = data.frame((Protein = data$Protein))
      colnames(matrix_final)[1] <- 'Protein'
      break
  }
}



for(i in 1:(length(name_data)/2)){
  data <- read.csv("data/MSstats_output.csv")
  
  temp.name = name_data[,i]
  data <- data[which(data$Label %in% temp.name), ]
  
  if (any(is.na(data))) {
      num = 1
      
      name_p = paste(name_data[,i+length(name_data)/2], '.p-value', sep = '')
      name_log2fc = paste(name_data[,i+length(name_data)/2], '.log2foldchange', sep = '')
      
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

write.csv(matrix_final, file = 'data/NAME_analytics.csv', quote = F, row.names = F)

}


