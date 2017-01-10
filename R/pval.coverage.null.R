#' Function to estimate the coverage of the overall change impact coefficient.
#'
#' @param newdat.ind simulated independent response data
#' @param model glm or gamMRSea model object
#' @param nsim number of simulations in newdat.ind
#' @param powercoefid coefficient id that relates to the eventchange term
#' @param empdistnull if known, provide the empirical distribution for the runs test.
#'
#' @author Lindesay Scott-Hayward
#'
#' @export

pval.coverage.null<-function(newdat.ind, newdat.corr=NULL, model, nsim, powercoefid, empdistnull=NULL){

  # if(!is.null(splineParams)){
  #   dists<-splineParams[[1]]$dist
  # }

  require(Matrix)
  require(mvtnorm)

  data<-model$data
  nc<-nsim/2
  nulldata.ind<-rbind(newdat.ind[data$eventphase==0, (1:nc)], newdat.ind[data$eventphase==0, ((nc+1):nsim)])
  if(is.null(newdat.corr)){
    nulldata.corr<-nulldata.ind
  }else{
    nulldata.corr<-rbind(newdat.corr[data$eventphase==0, (1:nc)], newdat.corr[data$eventphase==0, ((nc+1):nsim)])
  }
  nsim<-ncol(nulldata.ind)

  splineParams<-model$splineParams

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~ generate empirical distribution ~
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # calculate empirical distribution using the independent data
  if(is.null(empdistnull)){
    empdistnull<-getEmpDistribution(n.sim = nsim, simData=nulldata.ind,
                                      model = model, data = data, plot=FALSE,
                                      returnDist = TRUE, dots=FALSE)
  }

  imppvals= rawrob=vector(length=nsim)
  betacis<-matrix(NA, nrow=nsim, 2)

  for(i in 1:nsim){

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ fit model ~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # fit models using the correlated data (if available)
    data$response<-nulldata.corr[,i]
    sim_glm<-update(model, response ~ ., data=data)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~ power p-value ~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get runs test result using empirical distribution
    runspvalemp<-runsTest(residuals(sim_glm, type='pearson'),
                           emp.distribution = empdistnull)$p.value
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
      imppval<-summary(sim_glm)$rawp[powercoefid]
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

  } # end nsim
  return(list(imppvals=imppvals, betacis=betacis, empdistnull=empdistnull))
} # end function
