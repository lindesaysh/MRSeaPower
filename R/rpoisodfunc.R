rpois.od<-function(n, lambda, d=1){
  if(d[1]==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}
