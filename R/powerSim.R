powerSimOverallChange<-function(newdat, model, empdistribution, nsim, powercoefid, predictionGrid=NULL, g2k=NULL, splineParams=NULL, bootstrapCI=TRUE, sigdif=TRUE, n.boot=1000, impact.loc=NULL){

  require(mvtnorm)
  data<-model$data
  if(is.null(nsim)){nsim=ncol(newdat)}

  imppvals= rawrob=vector(length=nsim)
  powsimfits=matrix(NA, nrow=nrow(data), ncol=nsim)
  betacis<-matrix(NA, nrow=nsim, 2)
  #cis<-array(NA, c(nrow(predictionGrid), 2, nsim))
  indpvals = familypvals = list(nsim)
  preds<-matrix(NA, nrow=nrow(predictionGrid), ncol=nsim)
    bootdifferences=bootdifferences.pct=array(NA, c((nrow(predictionGrid)/2), n.boot, nsim))

  if(!is.null(impact.loc)){
    # make difference dataset
    preddiffs<-predictionGrid[predictionGrid$eventphase==0,]
    # distances from wfarm to each grid loc.
    preddiffs$impdist<-as.matrix(dist(x = rbind(impact.loc, preddiffs[,c('x.pos', 'y.pos')])))[1,-1]
    # cut the distances into bins
    br<-seq(0, max(preddiffs$impdist), length=10)
    preddiffs$cuts<-cut(preddiffs$impdist, breaks=br)
  }else{cat('No impact location given so distance to impact data not calculated.')}

  preddifferences<-matrix(NA, nrow=(nrow(predictionGrid)/2), ncol=nsim)

  for(i in 1:nsim){

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ fit model ~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    #dists<<-splineParams[[1]]$dists
    data$response<-newdat[,i]
    sim_glm<-update(model, response ~ ., data=data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ power p-value ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get runs test result using empirical distribution
    runspvalemp<-runs.test(residuals(sim_glm, type='pearson'),
                           critvals = empdistribution)$p.value

    sim.anv<-anova.gamMRSea(sim_glm)
    if(length(pmatch("LocalRadialFunction(radiusIndices, dists, radii, aR):eventphase", rownames(sim.anv)))>0){
      if(runspvalemp<=0.05){
        # significant, therefore correlated, use robust se
        class(sim_glm)<-c('gamMRSea', class(sim_glm))
        sim_glm$panels<-model$panels
        imppval<-sim.anv$P[nrow(sim.anv)]
        rawrob[i]<-1
        robust<-TRUE
      }else{
        # not significant, use raw se
        class(sim_glm)<-c('gamMRSea', class(sim_glm))
        sim_glm$panels<-1:nrow(data)
        imppval<-sim.anv$P[nrow(sim.anv)]
        rawrob[i]<-0
        robust<-FALSE
      }
    }else{
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
        class(sim_glm)<-c('gamMRSea', class(sim_glm))
        imppval<-summary(sim_glm)$rawp[powercoefid]
        rawrob[i]<-0
        robust<-FALSE
      }

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
    #dists<<-g2k
    preds[,i]<-predict.gamMRSea(predictionGrid, splineParams = splineParams, g2k = g2k, model=sim_glm, type = 'response')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ bootstrap CIs ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # if(bootstrapCI==T){
    #   bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid, splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)
    #     # # get upper and lower cis
    #     # quants<-c(0.025, 0.975)
    #     # cis[,,i]<-t(apply(bootPreds, 1,  quantile, probs= quants, na.rm=T ))
    #   } # end bootstrap

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Differences ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(sigdif==T){
      bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid, splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)

      # extract distribution of null differences
      nulldifferences<-bootPreds[predictionGrid$eventphase==0,(1:(n.boot/2))] - bootPreds[predictionGrid$eventphase==0,(((n.boot/2)+1):n.boot)]

      # find the predicted differences from this simulation
      preddifferences[,i]<-preds[predictionGrid$eventphase==1,i] - preds[predictionGrid$eventphase==0,i]

      # calculate p-values for individual and family wide cells.
      indpvals[[i]]<-pval.differences(nulldifferences, preddifferences[,i], family=FALSE)
      familypvals[[i]]<-pval.differences(nulldifferences, preddifferences[,i], family=TRUE)

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~ Abundance/mean proportion ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(sim_glm$family[[1]]=='poisson' | sim_glm$family[[1]]=='quasipoisson'){
        if(i==1){bsum=asum=vector(length=nsim)}
        bsum[i]<-sum(bootPreds[predictionGrid$eventphase==0,])
        asum[i]<-sum(bootPreds[predictionGrid$eventphase==1,])
      }
      if(sim_glm$family[[1]]=='binomial' | sim_glm$family[[1]]=='quasibinomial'){
        if(i==1){bmean=amean=vector(length=nsim)}
        bmean[i]<-mean(bootPreds[predictionGrid$eventphase==0,])
        amean[i]<-mean(bootPreds[predictionGrid$eventphase==1,])
      }

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~ Distance to Windfarm ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if(!is.null(impact.loc)){
        #if(nrow(predictionGrid)>400){
        preddifferences.pct<-(as.vector(preddifferences[,i])/as.vector(preds[predictionGrid$eventphase==0,i]))*100

        d2imp.plotdata<-matrix(NA, nrow=length(unique(preddiffs$cuts)), ncol=9)
        d2imp.plotdata[,1]<-tapply(preddiffs$impdist, preddiffs$cuts, mean)
        d2imp.plotdata[,2]<-tapply(preddifferences[,i], preddiffs$cuts, mean)
        d2imp.plotdata[,6]<-tapply(preddifferences.pct, preddiffs$cuts, mean)
        colnames(d2imp.plotdata)<-c('MeanDist', 'MeanDiff', 'bootMeanDiff', 'LowerCI', 'UpperCI', 'MeanDiff.pct', 'bootMeanDiff.pct', 'LowerCI.pct', 'UpperCI.pct')

        bootdifferences[,,i]<-bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,]
        bootdifferences.pct[,,i]<-((bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,])/bootPreds[predictionGrid$eventphase==0,])*100
      }

    } # end sigdif

  } # end nsim

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Abundance/mean proportion ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(sim_glm$family[[1]]=='poisson' | sim_glm$family[[1]]=='quasipoisson'){
      quants<-c(0.025, 0.975)
      abund<-matrix(NA, 2, 3)
      abund[,1]<-c(mean(bsum), mean(asum))
      abund[1,2:3]<-quantile(bsum, probs=quants)
      abund[2,2:3]<-quantile(asum, probs=quants)
      rownames(abund)<-c('Before', 'After')
      colnames(abund)<-c('Abundance', 'LowerCI', 'UpperCI')
    }
    if(sim_glm$family[[1]]=='binomial' | sim_glm$family[[1]]=='quasibinomial'){
      quants<-c(0.025, 0.975)
      meanp<-matrix(NA, 2, 3)
      meanp[,1]<-c(mean(bmean), mean(amean))
      meanp[1,2:3]<-quantile(bmean, probs=quants)
      meanp[2,2:3]<-quantile(amean, probs=quants)
      rownames(meanp)<-c('Before', 'After')
      colnames(meanp)<-c('Mean', 'LowerCI', 'UpperCI')
    }


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Distance to Windfarm ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(!is.null(impact.loc)){

      # get the upper and lower quantile for the differences in each cut point.
      uniquecuts<-unique(preddiffs$cuts)
      for(p in 1:length(uniquecuts)){
        d2imp.plotdata[p,4:5]<-quantile(na.omit(bootdifferences[preddiffs$cuts==uniquecuts[p],,]), probs=c(0.025, 0.975))
        d2imp.plotdata[p,3]<-mean(na.omit(bootdifferences[preddiffs$cuts==uniquecuts[p],,]))
        d2imp.plotdata[p,8:9]<-quantile(na.omit(bootdifferences.pct[preddiffs$cuts==uniquecuts[p],,]), probs=c(0.025, 0.975))
        d2imp.plotdata[p,7]<-mean(na.omit(bootdifferences.pct[preddiffs$cuts==uniquecuts[p],,]))
      }

    } # end impact.loc



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Return list object~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  output<-list(rawrob=rawrob, imppvals=imppvals, betacis=betacis, powsimfits=powsimfits, significant.differences=list(individual=indpvals, family=familypvals), d2imp.plotdata=d2imp.plotdata)

  if(sim_glm$family[[1]]=='poisson' | sim_glm$family[[1]]=='quasipoisson'){
    output$Abundance = abund
  }
  if(sim_glm$family[[1]]=='binomial' | sim_glm$family[[1]]=='quasibinomial'){
    output$Mean.proportion = meanp
  }

  class(output)<-'gamMRSea.power'

  return(output)
}

