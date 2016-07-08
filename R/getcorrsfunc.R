acffunc_dat<-function (block, data) 
{
  blocktab <- table(block)
  acfmat <- matrix(NA, length(unique(block)), max(blocktab))
  for (i in 1:length(unique(block))) {
    print(i)
    corr <- as.vector(acf(data[which(block == unique(block)[i])], plot = F, lag.max = max(blocktab))$acf)
    if(length(which(is.na(corr)))>0){
      corr<-c(1, rep(0.999999, (length(corr)-1)))
      #corr<-seq(1, 0.001, length.out=(length(corr)))
    }
    acfmat[i, 1:length(corr)] <- corr
  }
  return(list(acfmat = acfmat, blocktab = blocktab))
}