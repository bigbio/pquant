getSelector <- function(fileData, flag, DDA2009.proposed){
      # Automatically create the manually created matrix in MSstats, user manual p23
      len <- length(levels(DDA2009.proposed$FeatureLevelData$GROUP))
      
      tmp <- t(combn(len,2))
      matrix_len = length(t(combn(len,2))) / 2
      
      ourMatrix <- matrix(c(0:0),nrow=matrix_len,ncol=len)
      
      for(i in 1:matrix_len){
        ourMatrix[i, tmp[i]] = -1
        ourMatrix[i, tmp[i + matrix_len]] = 1
      }
      
      ourCondition <- levels(DDA2009.proposed$ProteinLevelData$GROUP)
      tmp_name <- matrix(ourCondition, nr=len, nc=1)
      name <- matrix(nr=matrix_len, nc=1)
      for(i in 1:matrix_len){
        name[i,1] <- sprintf('%s-%s', tmp_name[tmp[i+matrix_len]], tmp_name[tmp[i]])
      }
      
      selector <- name
      
      return(selector)

}