powerSimRunsTest<-function(newdat, model, nsim, powercoefid){
  
  require(mvtnorm)
  data<-model$data
  if(is.null(nsim)){nsim=ncol(newdat)}
  
  imppvals=matrix(NA, nsim, 4) 
  rawrob=vector(length=nsim)
  betacis<-matrix(NA, nrow=nsim, 2)
  
  
  for(i in 1:nsim){
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ fit model ~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    #dists<<-splineParams[[1]]$dists
    data$response<-newdat[,i]
    sim_glm<-update(model, response ~ ., data=data)
    sim_gee<-geeglm(sim_glm$formula, data=data, family='poisson', id = panels, corstr = 'ar1')
    sim_geeind<-geeglm(sim_glm$formula, data=data, family='poisson', id = panels)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ power p-value ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(sim_glm)<-c('gamMRSea', class(sim_glm))
    sim_glm$panels<-data$panels
    imppval.rob<-summary(sim_glm)$coefficients[powercoefid,5]
    imppval.raw<-summary(sim_glm)$rawp[powercoefid]
    imppval.gee<-summary(sim_gee)$coefficients[powercoefid,4]
    imppval.geeind<-summary(sim_geeind)$coefficients[powercoefid,4]
    
    imppvals[i,]<-c(imppval.raw, imppval.rob, imppval.gee, imppval.geeind)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~ beta CI's ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    #   est<- coefficients(sim_glm)
    #   if(robust==T){
    #     vbeta<-summary(sim_glm)$cov.robust
    #     samplecoeff<-NULL
    #     try(samplecoeff<- rmvnorm(500,est,vbeta, method='svd'))
    #     if(is.null(samplecoeff)){
    #       vbeta<-as.matrix(nearPD(as.matrix(summary(sim_glm)$cov.robust))$mat)
    #       samplecoeff<- rmvnorm(500,est,vbeta, method='svd')
    #     }
    #   }
    #   else{
    #     vbeta<-summary(sim_glm)$cov.scaled
    #     samplecoeff<-NULL
    #     try(samplecoeff<- rmvnorm(500,est,vbeta, method='svd'))
    #     if(is.null(samplecoeff)){
    #       vbeta<-as.matrix(nearPD(as.matrix(summary(sim_glm)$cov.scaled))$mat)
    #       samplecoeff<- rmvnorm(500,est,vbeta, method='svd')
    #     }
    #   }
    #   betacis[i,]<- as.numeric(quantile(samplecoeff[,powercoefid],probs = c(0.025, 0.975)))
  } 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Return list object~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  output<-list(rawrob=rawrob, imppvals=imppvals, betacis=betacis)
  
  class(output)<-'gamMRSea.power'
  
  return(output)
}
