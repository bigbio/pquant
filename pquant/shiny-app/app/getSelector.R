getSelector <- function(fileData, flag, DDA2009.proposed, DDA2009.TMP){

  if (flag == 'qc'){
      tmp <- levels(DDA2009.proposed$ProcessedData$PROTEIN)
      selector <- append('allonly', tmp, 1)
  }
  
  else {
      # Automatically create the manually created matrix in MSstats, user manual p23
      len <- length(levels(DDA2009.TMP$FeatureLevelData$GROUP))
      
      ourMatrix <- matrix(c(0:0),nrow=len,ncol=len)
      diag(ourMatrix) = -1
      for(i in 1:len-1){
        ourMatrix[i,i+1] = 1
      }
      ourMatrix[len,1] = 1
      
      ourCondition <- levels(DDA2009.TMP$ProteinLevelData$GROUP)
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
