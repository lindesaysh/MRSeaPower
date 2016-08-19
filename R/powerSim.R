powerSimOverallChange<-function(newdat, model, empdistribution, nsim, powercoefid, predictionGrid=NULL, g2k=NULL, splineParams=NULL, bootstrapCI=TRUE, sigdif=TRUE, n.boot=1000){

  require(mvtnorm)
  data<-model$data
  if(is.null(nsim)){nsim=ncol(newdat)}

  imppvals= rawrob=vector(length=nsim)
  powsimfits=matrix(NA, nrow=nrow(data), ncol=nsim)
  betacis<-matrix(NA, nrow=nsim, 2)
  #cis<-array(NA, c(nrow(predictionGrid), 2, nsim))
  indpvals = familypvals = list(nsim)
  preds<-matrix(NA, nrow=nrow(predictionGrid), ncol=nsim)

  for(i in 1:nsim){

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ fit model ~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    data$response<-newdat[,i]
    sim_glm<-update(model, response ~ ., data=data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ power p-value ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get runs test result using empirical distribution
    runspvalemp<-runs.test(residuals(sim_glm, type='pearson'),
                           critvals = empdistribution)$p.value
    if(runspvalemp<=0.05){
      # significant, therefore correlated, use robust se
      class(sim_glm)<-c('gamMRSea', class(sim_glm))
      sim_glm$panels<-data$panelid
      imppval<-summary(sim_glm)$coefficients[powercoefid,5]
      #class(sim_glm)<-class(sim_glm)[-1]
      rawrob[i]<-1
      robust<-TRUE
    }else{
      # not significant, use raw se
      imppval<-summary(sim_glm)$coefficients[powercoefid,4]
      rawrob[i]<-0
      robust<-FALSE
    }
    imppvals[i]<-imppval

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~ beta CI's ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    est<- coefficients(sim_glm)
    if(robust==T){
      vbeta<-summary(sim_glm)$cov.robust
      samplecoeff<-NULL
      try(samplecoeff<- rmvnorm(500,est,vbeta, method='svd'))
      if(is.null(samplecoeff)){
        vbeta<-as.matrix(nearPD(as.matrix(summary(sim_glm)$cov.robust))$mat)
        samplecoeff<- rmvnorm(500,est,vbeta, method='svd')
      }
    }
    else{
      vbeta<-summary(sim_glm)$cov.scaled
      samplecoeff<-NULL
      try(samplecoeff<- rmvnorm(500,est,vbeta, method='svd'))
      if(is.null(samplecoeff)){
        vbeta<-as.matrix(nearPD(as.matrix(summary(sim_glm)$cov.scaled))$mat)
        samplecoeff<- rmvnorm(500,est,vbeta, method='svd')
      }
    }
    betacis[i,]<- as.numeric(quantile(samplecoeff[,powercoefid],probs = c(0.025, 0.975)))
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ fitted values ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    powsimfits[,i]<-fitted(sim_glm)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ predictions to grid ~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    preds[,i]<-predict.gamMRSea(predictionGrid, splineParams = splineParams, g2k = g2k, model=sim_glm, type = 'response')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ bootstrap CIs ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # if(bootstrapCI==T){
    #   bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid, splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)
    #
    #     # # get upper and lower cis
    #     # quants<-c(0.025, 0.975)
    #     # cis[,,i]<-t(apply(bootPreds, 1,  quantile, probs= quants, na.rm=T ))
    #
    #   } # end bootstrap
    #

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Differences ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(sigdif==T){
      bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid[predictionGrid$eventphase==0,], splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)

      nulldifferences<-bootPreds[,(1:(n.boot/2))] - bootPreds[,(((n.boot/2)+1):n.boot)]
      preddifferences<-preds[predictionGrid$eventphase==1,i] - preds[predictionGrid$eventphase==0,i]
      indpvals[[i]]<-pval.differences(nulldifferences, preddifferences, family=FALSE)
      familypvals[[i]]<-pval.differences(nulldifferences, preddifferences, family=TRUE)
    # # get upper and lower cis
    # quants<-c(0.025, 0.975)
    # cis[,,i]<-t(apply(bootPreds, 1,  quantile, probs= quants, na.rm=T ))


    } # end sigdif

  } # end nsim


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Return list object~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  return(list(rawrob=rawrob, imppvals=imppvals, betacis=betacis, powsimfits=powsimfits, significant.differences=list(inividual=indpvals, family=familypvals)))
}

