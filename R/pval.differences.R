#' Function to get p-values for individual cell differences
#'
#'
#' @param null
#' @param estimate
#'
#' @example
#'

pval.differences<-function(null, estimate, family=FALSE){

  if(family==FALSE){
    if(length(estimate)>1){
      if(nrow(null)!=length(estimate)){warning('Null distribution and estimate dimensions do not match')}
      onep=twop=vector(length=length(estimate))
      for(i in 1:length(estimate)){
        onep[i]<-nulldistpval(null[i,], estimate[i], type='one')
        twop[i]<-nulldistpval(null[i,], estimate[i], type='two')
      }
    }else{
      onep<-nulldistpval(null, estimate, type='one')
      twop<-nulldistpval(null, estimate, type='two')
    }
  }
  if(family==TRUE){
    if(length(estimate)>1){
      if(nrow(null)!=length(estimate)){warning('Null distribution and estimate dimensions do not match')}
      onep=twop=vector(length=length(estimate))
      for(i in 1:length(estimate)){
        onep[i]<-nulldistpval(null, estimate[i], type='one')
        twop[i]<-nulldistpval(null, estimate[i], type='two')
      }
    }else{
      onep<-nulldistpval(null, estimate, type='one')
      twop<-nulldistpval(null, estimate, type='two')
    }
  }

  return(list(onesided = onep, twosided = twop))

}


nulldistpval<-function(distribution, est, type='two'){

  if(type=='one'){
    if(est<mean(distribution)){
      l<-length(which(distribution<=est))
      p<-l/length(distribution)
    }else{
      l<-length(which(distribution>=est))
      p<-l/length(distribution)
    }
  }
  if(type=='two'){
    l<-length(which(distribution<=-abs(est))) + length(which(distribution>=abs(est)))
    p<-l/length(distribution)
  }
  return(p)
}


