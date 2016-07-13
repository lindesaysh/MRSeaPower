
getRunsCritVals<-function(n.sim, simData, model, data, plot=FALSE, returnDist=FALSE){
  runsstatH0<-vector(length=n.sim)
for(i in 1:n.sim){
  if((i/500)%%1 == 0){cat(i, '\n')}else{cat('.')}
  data$response=simData[,i]
  sim_glm<- update(model, response ~ ., data=data)
  runs<-runs.test(residuals(sim_glm, type='pearson'))
  runsstatH0[i]<-runs$statistic
}

  critvals<-quantile(runsstatH0, probs = c(0.025, 0.975))

  if(plot){
    hist(runsstatH0, main='Empirical Distribution: Runs Test Statistic', xlim=range(c(qnorm(0.025), qnorm(0.975), runsstatH0)), xlab='Test Statistic')
    abline(v=critvals, col='red', lwd=2)
    abline(v=c(qnorm(0.025), qnorm(0.975)), col='blue', lwd=2)
  }

  if(returnDist){output<-runsstatH0}else{output<-as.vector(critvals)}

  return(output)
}
