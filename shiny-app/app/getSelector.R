getSelector <- function(fileData, flag){
  #default option
  DDA2009.proposed <- dataProcess(raw = fileData,
                                  normalization = 'equalizeMedians',
                                  summaryMethod = 'TMP',
                                  censoredInt = "NA",
                                  cutoffCensored = "minFeature",
                                  MBimpute = TRUE,
                                  maxQuantileforCensored=0.999)
  if (flag == 'qc'){
      tmp <- levels(DDA2009.proposed$ProcessedData$PROTEIN)
      selector <- append('allonly', tmp, 1)
  }
  
  else if (flag == 'volcano'){
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
      
      selector <- name
  }
}