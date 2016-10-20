
getRunsCritVals.raw<-function(n.sim, simData, model, data, plot=FALSE, returnDist=FALSE,dots=TRUE){
  runsstatH0=betas=vector(length=n.sim)
  model$panels<-1:nrow(data) # make sure that independent panels specified for null distribution
  for(i in 1:n.sim){
    if(dots==TRUE){if((i/500)%%1 == 0){cat(i, '\n')}else{cat('.')}}
    
    splineParams = model$splineParams
    data$response=simData[,i]
    sim_glm<- update(model, response ~ ., data=data)
    
    d<-summary(sim_glm)$dispersion
    #print(d)
    resid<-residuals(sim_glm, type="response")
    
    runs<-runs.test(resid)
    runsstatH0[i]<-runs$statistic
    # if(imp.beta.coverage){
    #   betas[i]<-coefficients(model)[length(coefficients(model))]
    # }
  }
  
  teststat<-na.omit(runsstatH0)
  if(any(abs(teststat)=='Inf')){
    teststat<-teststat[-which(abs(teststat) == 'Inf')]
  }
  
  critvals<-quantile(teststat, probs = c(0.025, 0.975))
  
  if(plot){
    hist(teststat, main='Empirical Distribution: Runs Test Statistic', xlim=range(c(qnorm(0.025), qnorm(0.975), teststat)), xlab='Test Statistic')
    abline(v=critvals, col='red', lwd=2)
    abline(v=c(qnorm(0.025), qnorm(0.975)), col='blue', lwd=2)
  }
  
  
  
  # if(returnDist){
  #   output<-list=c(empirical.distribution.null=runsstatH0, beta.coverage.null = )
  # }else{
  #   output<-list=c(empirical.critical=as.vector(critvals), beta.coverage.null = )
  #   }
  if(returnDist){
    output<-runsstatH0
  }else{
    output<-as.vector(critvals)
  }
  return(output)
}
