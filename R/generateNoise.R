
generateNoise<-function(n, response, family, ...){
  simData<-matrix(NA, nrow=length(response), ncol=n)
  for(i in 1:n){
    if(family=='poisson'){
      simData[,i]<-rpois.od(length(response), response, ...)    
    }
    if(family=='binomial'){
      simData[,i]<-rbinom(n=length(response), prob=response, ...)
    }
  }
  return(simData)
}