runsPowerSim<-function(newdat, model, empdistribution, nsim, powercoefid){

  data<-model$data
  if(is.null(nsim)){nsim=ncol(newdat)}

  imppvals= rawrob=matrix(NA, nsim, 2)
  for(i in 1:nsim){
    data$response<-newdat[,i]
    sim_glm<-update(model, response ~ ., data=data)

    # get runs test result using empirical distribution
    runspvalemp<-runs.test(residuals(sim_glm, type='pearson'),
                           critvals = empdistribution)$p.value
    if(runspvalemp<=0.05){

      # significant, therefore correlated, use robust se
      class(sim_glm)<-c('gamMRSea', class(sim_glm))
      sim_glm$panels<-data$panelid
      imppval<-summary(sim_glm)$coefficients[powercoefid,5]
      class(sim_glm)<-class(sim_glm)[-1]
      rawrob[i,1]<-1
    }else{
      # not significant, use raw se
      imppval<-summary(sim_glm)$coefficients[powercoefid,4]
      rawrob[i,1]<-0
    }
    imppvals[i,1]<-imppval

    # traditional runs test
    runspvalN<-runs.test(residuals(sim_glm, type='pearson'))$p.value
    if(runspvalN<=0.05){
      #print('N: Not correlated')
      # significant, robust se
      class(sim_glm)<-c('gamMRSea', class(sim_glm))
      sim_glm$panels<-data$panelid
      imppvalN<-summary(sim_glm)$coefficients[powercoefid,5]
      rawrob[i,2]<-1
    }else{
      # not significant, raw se
      imppvalN<-summary(sim_glm)$coefficients[powercoefid,4]
      rawrob[i,2]<-0
    }
    imppvals[i,2]<-imppvalN
  }
  colnames(rawrob)= colnames(imppvals) = c('Emp', 'Norm')
  return(list(rawrob=rawrob, imppvals=imppvals))
}

