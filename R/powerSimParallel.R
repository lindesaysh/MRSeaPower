#' Function to run a power analysis for a model of type 'gamMRSea'
#'
#' @param newdat matrix of values whose rows are observations and columns are simulated sets of response data
#' @param model model object of class `gamMRSea`.  This is the model used to generate the data, updated for the inclusion of an event term.
#' @param empdistribution Empirical distribution for the runs test
#' @param nsim number of simulations to compute
#' @param powercoefid the ID of the coefficient in \code{model} that represents the event term
#' @param predictionGrid dataframe containing a coordinate grid with covariates.  This grid is used in the assessment of before and after change differences.
#' @param g2k object required for \code{MRSea} model with spline terms
#' @param splineParams object required for \code{MRSea} model with spline terms
#' @param n.boot Number of bootstraps to use to get percentile based confidence intervals for predictions.
#' @param nCores How many cores to use for parallel processing.  If one, then parallel computing is not used and a progress bar is printed. If parallel, estimated time till finish along with a progress bar is provided.
#'
#' @return An object of class \code{gamMRSea.power}.
#'
#' \item{rawrob}{Vector of zeros and ones indicating whether raw or robust s.e. have been used.}
#' \item{imppvals}{p-values for the event change term for each simulated run}
#' \item{betacis}{Confidence intervals for the event term coefficient (only if site-wide change)}
#' \item{powsimfits}{Fitted values each model}
#' \item{bootdifferences}{data frame of bootstrapped differences (before and after event). Number of rows is the number of rows in the prediction grid divided by 2}
#' \item{bootpreds}{data frame of bootstrapped predictions.}
#' \item{Abundance}{Table of sitewide abundance (with upper and lower 95\% CI)}
#'
#'
#' @examples
#' # See the vignette for an example of how to run this code.
#'
#' @author LAS Scott-Hayward, University of St Andrews
#'
#' @export
#'
#'
powerSimPll<-function(newdat, model, empdistribution, nsim, powercoefid, predictionGrid=NULL, g2k=NULL, splineParams=NULL, n.boot=500, nCores=1){

  sigdif<-TRUE

  require(mvtnorm)
  data<-model$data
  if(is.null(nsim)){nsim=ncol(newdat)}

  imppvals= rawrob=vector(length=nsim)
  powsimfits=matrix(NA, nrow=nrow(data), ncol=nsim)
  betacis<-matrix(NA, nrow=nsim, 2)
  #cis<-array(NA, c(nrow(predictionGrid), 2, nsim))
  indpvals = familypvals = list(nsim)
  preds<-matrix(NA, nrow=nrow(predictionGrid), ncol=nsim)
  bootdifferences=array(NA, c((nrow(predictionGrid)/2), n.boot, nsim))

  # if(!is.null(impact.loc)){
  #   # make difference dataset
  #   preddiffs<-predictionGrid[predictionGrid$eventphase==0,]
  #   # distances from wfarm to each grid loc.
  #   preddiffs$impdist<-as.matrix(dist(x = rbind(impact.loc, preddiffs[,c('x.pos', 'y.pos')])))[1,-1]
  #   # cut the distances into bins
  #   br<-seq(0, max(preddiffs$impdist), length=10)
  #   preddiffs$cuts<-cut(preddiffs$impdist, breaks=br)
  # }else{cat('No impact location given so distance to impact data not calculated.')}

  preddifferences<-matrix(NA, nrow=(nrow(predictionGrid)/2), ncol=nsim)



  if(nCores>1){
    require(parallel)
    cat('Code in parallel')

    computerCores <- getOption("cl.cores", detectCores())
    if(nCores>(computerCores-2)){nCores<-(computerCores-2)}
    myCluster <- makeCluster(nCores) ; myCluster # initialise the cluster
    clusterExport(myCluster, ls(), envir=environment()) # export all objects to each cluster
    # export directory and functions to each cluster
    clusterEvalQ(myCluster, {
      require(splines)
      require(mgcv)
      require(MRSea)
      require(dplyr)
      require(Matrix)
      require(mvtnorm)
      require(MRSeaPower)
    })

    # only do parametric boostrap if no data re-sampling and no nhats provided
    Routputs<-pbapply::pblapply(cl=myCluster, X=1:nsim, FUN=function(i){

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
      runspvalemp<-runsTest(residuals(sim_glm, type='pearson'),
                             emp.distribution = empdistribution)$p.value

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

      imppvals<-imppval

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
      betacis<- as.numeric(quantile(samplecoeff[,powercoefid],probs = c(0.025, 0.975)))
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~ fitted values ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      powsimfits<-fitted(sim_glm)

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~ predictions to grid ~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      #dists<<-g2k
      preds<-predict.gamMRSea(predictionGrid, g2k = g2k, model=sim_glm, type = 'response')

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~ Differences ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(sigdif==T){
        bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid, splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)


        # find the predicted differences from this simulation
        preddifferences<-preds[predictionGrid$eventphase==1] - preds[predictionGrid$eventphase==0]

        bootdifferences<-bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,]
        #bootdifferences.pct<-((bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,])/bootPreds[predictionGrid$eventphase==0,])*100


        # # calculate p-values for individual and family wide cells.
        # indpvals<-pval.differences(nulldifferences, preddifferences, family=FALSE)
        # familypvals<-pval.differences(nulldifferences, preddifferences, family=TRUE)
        #
        # indpvals<-as.matrix(data.frame(indpvals))
        # familypvals<-as.matrix(data.frame(familypvals))

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~ Abundance/mean proportion ~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # remove bootstraps with excessive predictions.
        badids<-which(colMeans(bootPreds)>max(preds))

        if(length(badids)>0){
          for(bad in 1:length(badids)){
            bootPreds[, badids[bad]]<-NA
            bootdifferences[, badids[bad]]<-NA
          }
        }

        if(sim_glm$family[[1]]=='poisson' | sim_glm$family[[1]]=='quasipoisson'){
          #if(i==1){bsum=asum=vector(length=nsim)}
          bsum<-sum(bootPreds[predictionGrid$eventphase==0,])
          asum<-sum(bootPreds[predictionGrid$eventphase==1,])
        }
        if(sim_glm$family[[1]]=='binomial' | sim_glm$family[[1]]=='quasibinomial'){
          #if(i==1){bmean=amean=vector(length=nsim)}
          bsum<-mean(bootPreds[predictionGrid$eventphase==0,])
          asum<-mean(bootPreds[predictionGrid$eventphase==1,])
        }

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~ Distance to Windfarm ~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # if(!is.null(impact.loc)){
        #   #if(nrow(predictionGrid)>400){
        #   preddifferences.pct<-(as.vector(preddifferences)/as.vector(preds[predictionGrid$eventphase==0]))*100
        #
        #   d2imp.plotdata<-matrix(NA, nrow=length(unique(preddiffs$cuts)), ncol=9)
        #   d2imp.plotdata[,1]<-tapply(preddiffs$impdist, preddiffs$cuts, mean)
        #   d2imp.plotdata[,2]<-tapply(preddifferences, preddiffs$cuts, mean)
        #   d2imp.plotdata[,6]<-tapply(preddifferences.pct, preddiffs$cuts, mean)
        #   colnames(d2imp.plotdata)<-c('MeanDist', 'MeanDiff', 'bootMeanDiff', 'LowerCI', 'UpperCI', 'MeanDiff.pct', 'bootMeanDiff.pct', 'LowerCI.pct', 'UpperCI.pct')
        #

        # }

      } # end sigdif


      return(list(imppvals=imppvals, betacis=betacis, powsimfits=powsimfits, preds=preds, bsum=bsum, asum=asum, bootPreds=bootPreds, bootdifferences=bootdifferences, preddifferences=preddifferences))
    })

    stopCluster(myCluster)

    cat("Creating Outputs...")

    imppvals= sapply(Routputs, '[[','imppvals')
    bsum= sapply(Routputs, '[[','bsum')
    asum= sapply(Routputs, '[[','asum')
    powsimfits=sapply(Routputs, '[[','powsimfits')
    betacis=t(sapply(Routputs, '[[','betacis'))
    # indpvals = sapply(Routputs, '[[','indpvals', simplify='array')
    # familypvals = sapply(Routputs, '[[','familypvals', simplify='array')
    preds<-sapply(Routputs, '[[','preds')
    preddifferences<-sapply(Routputs, '[[','preddifferences')
    bootdifferences=sapply(Routputs, '[[','bootdifferences', simplify='array')
    #bootdifferences.pct=sapply(Routputs, '[[','bootdifferences.pct', simplify  ='array')
    #d2imp.plotdata=Routputs[[nsim]]$d2imp.plotdata
    bootPreds.all=sapply(Routputs, '[[','bootPreds', simplify='array')

    rm(Routputs)
    gc()

   #~~~
    # try alply
    # require(plyr)
    # indpvals<-alply(indpvals,3)
    # familypvals<-alply(familypvals,3)
    # detach(package:plyr)
    #~~~

    }else{

      bootPreds.all<-array(NA, c(nrow(predictionGrid), n.boot,  nsim))
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
    runspvalemp<-runsTest(residuals(sim_glm, type='pearson'),
                           emp.distribution = empdistribution)$p.value

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
    preds[,i]<-predict.gamMRSea(predictionGrid, g2k = g2k, model=sim_glm, type = 'response')

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Differences ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(sigdif==T){
      bootPreds<-do.bootstrap.cress.robust(sim_glm,predictionGrid, splineParams=splineParams, g2k=g2k, B=n.boot, robust=robust, cat.message=FALSE)

      # # extract distribution of null differences
      # nulldifferences<-bootPreds[predictionGrid$eventphase==0,(1:(n.boot/2))] - bootPreds[predictionGrid$eventphase==0,(((n.boot/2)+1):n.boot)]

      # find the predicted differences from this simulation
      preddifferences[,i]<-preds[predictionGrid$eventphase==1,i] - preds[predictionGrid$eventphase==0,i]


      bootdifferences[,,i]<-bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,]
      #bootdifferences.pct[,,i]<-((bootPreds[predictionGrid$eventphase==1,] - bootPreds[predictionGrid$eventphase==0,])/bootPreds[predictionGrid$eventphase==0,])*100

      # # calculate p-values for individual and family wide cells.
      # indpvals[[i]]<-pval.differences(nulldifferences, preddifferences[,i], family=FALSE)
      # familypvals[[i]]<-pval.differences(nulldifferences, preddifferences[,i], family=TRUE)

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~ Abundance/mean proportion ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # remove bootstraps with excessive predictions.
      badids<-which(colMeans(bootPreds)>max(preds))

      if(length(badids)>0){
        for(bad in 1:length(badids)){
          bootPreds[, badids[bad]]<-NA
          bootdifferences[, badids[bad],i]<-NA
        }
      }

      if(sim_glm$family[[1]]=='poisson' | sim_glm$family[[1]]=='quasipoisson'){
        if(i==1){bsum=asum=vector(length=nsim)}
        bsum[i]<-sum(bootPreds[predictionGrid$eventphase==0,])
        asum[i]<-sum(bootPreds[predictionGrid$eventphase==1,])
      }
      if(sim_glm$family[[1]]=='binomial' | sim_glm$family[[1]]=='quasibinomial'){
        if(i==1){bsum=asum=vector(length=nsim)}
        bsum[i]<-mean(bootPreds[predictionGrid$eventphase==0,])
        asum[i]<-mean(bootPreds[predictionGrid$eventphase==1,])
      }

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # ~~~~ Distance to Windfarm ~~~~~~~~~
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      # if(!is.null(impact.loc)){
      #   #if(nrow(predictionGrid)>400){
      #   preddifferences.pct<-(as.vector(preddifferences[,i])/as.vector(preds[predictionGrid$eventphase==0,i]))*100
      #
      #   d2imp.plotdata<-matrix(NA, nrow=length(unique(preddiffs$cuts)), ncol=9)
      #   d2imp.plotdata[,1]<-tapply(preddiffs$impdist, preddiffs$cuts, mean)
      #   d2imp.plotdata[,2]<-tapply(preddifferences[,i], preddiffs$cuts, mean)
      #   d2imp.plotdata[,6]<-tapply(preddifferences.pct, preddiffs$cuts, mean)
      #   colnames(d2imp.plotdata)<-c('MeanDist', 'MeanDiff', 'bootMeanDiff', 'LowerCI', 'UpperCI', 'MeanDiff.pct', 'bootMeanDiff.pct', 'LowerCI.pct', 'UpperCI.pct')
      #

         bootPreds.all[,,i]<-bootPreds
      # }

    } # end sigdif

  } # end nsim
} # end cores else statement

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Abundance/mean proportion ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(model$family[[1]]=='poisson' | model$family[[1]]=='quasipoisson'){
      quants<-c(0.025, 0.975)
      abund<-matrix(NA, 2, 3)
      abund[,1]<-c(mean(na.omit(bsum)), mean(na.omit(asum)))
      abund[1,2:3]<-quantile(bsum, probs=quants, na.rm=T)
      abund[2,2:3]<-quantile(asum, probs=quants, na.rm=T)
      rownames(abund)<-c('Before', 'After')
      colnames(abund)<-c('Abundance', 'LowerCI', 'UpperCI')
    }
    if(model$family[[1]]=='binomial' | model$family[[1]]=='quasibinomial'){
      quants<-c(0.025, 0.975)
      meanp<-matrix(NA, 2, 3)
      meanp[,1]<-c(mean(na.omit(bsum)), mean(na.omit(asum)))
      meanp[1,2:3]<-quantile(bsum, probs=quants, na.rm=T)
      meanp[2,2:3]<-quantile(asum, probs=quants, na.rm=T)
      rownames(meanp)<-c('Before', 'After')
      colnames(meanp)<-c('Mean', 'LowerCI', 'UpperCI')
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ null differences and p-vals ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nulldifferences<-preds[predictionGrid$eventphase==0,(1:(nsim/2))] - preds[predictionGrid$eventphase==0,(((nsim/2)+1):nsim)]

    for(pval in 1:nsim){
      indpvals[[pval]]<-pval.differences(nulldifferences, preddifferences[, pval], family=FALSE)
      #familypvals[[pval]]<-pval.differences(nulldifferences, preddifferences[, pval], family=TRUE)
    }


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~ Distance to Windfarm ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # if(!is.null(impact.loc)){
    #
    #   # get the upper and lower quantile for the differences in each cut point.
    #   uniquecuts<-unique(preddiffs$cuts)
    #   for(p in 1:length(uniquecuts)){
    #     d2imp.plotdata[p,4:5]<-quantile(na.omit(bootdifferences[preddiffs$cuts==uniquecuts[p],,]), probs=c(0.025, 0.975))
    #     d2imp.plotdata[p,3]<-mean(na.omit(bootdifferences[preddiffs$cuts==uniquecuts[p],,]))
    #     d2imp.plotdata[p,8:9]<-quantile(na.omit(bootdifferences.pct[preddiffs$cuts==uniquecuts[p],,]), probs=c(0.025, 0.975))
    #     d2imp.plotdata[p,7]<-mean(na.omit(bootdifferences.pct[preddiffs$cuts==uniquecuts[p],,]))
    #   }

    #} # end impact.loc

# remove bad bootstraps

badids<-which(colMeans(bootPreds.all)>max(preds), arr.ind=T)

if(length(badids)>0){
  for(bad in 1:nrow(badids)){
    bootPreds.all[, badids[bad,1], badids[bad,2]]<-NA
    bootdifferences[, badids[bad,1], badids[bad,2]]<-NA
  }
}

bootpredsmean<-rowMeans(bootPreds.all, na.rm=TRUE)
bootpredscis<-t(apply(bootPreds.all, 1, quantile, probs=c(0.025, 0.975), na.rm=TRUE))
estpreds<-data.frame(mean=bootpredsmean, bootpredscis,predictionGrid[,c('x.pos', 'y.pos', 'eventphase')])

bootdiffmean<-rowMeans(bootdifferences, na.rm=TRUE)
bootdiffcis<-t(apply(bootdifferences, 1, quantile, probs=c(0.025, 0.975), na.rm=TRUE))
estdiffs<-data.frame(mean=bootdiffmean, bootdiffcis, predictionGrid[predictionGrid$eventphase==0,c('x.pos', 'y.pos')])



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ Return list object~~~~~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~
  output<-list(rawrob=rawrob, imppvals=imppvals, betacis=betacis, powsimfits=powsimfits, significant.differences=list(individual=indpvals), bootdifferences=estdiffs, bootpreds=estpreds)

  if(model$family[[1]]=='poisson' | model$family[[1]]=='quasipoisson'){
    output$Abundance = abund
  }
  if(model$family[[1]]=='binomial' | model$family[[1]]=='quasibinomial'){
    output$Mean.proportion = meanp
  }

  class(output)<-'gamMRSea.power'

  return(output)
}

