# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

generateDistribData<-function(nsim, model, coeff, data, dist.func=NULL){
  rcoef<-rmvnorm(nsim, coeff, summary(model)$cov.unscaled)
  newdata<-model$family$linkinv(model.matrix(model)%*%t(rcoef))
  
  # make new datasets from fitted values of gee or gam:
  newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
  for(i in 1:nsim){
    newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata[,i], d = as.numeric(summary(model)$dispersion[1])) 
  }
  return(newdata2ndsim)
}

generateDistribData2<-function(nsim, model, newdata, betasamp=FALSE, coeff=NULL, dist.func='quasipoisson', disp){
  
  if(betasamp==TRUE){
    rcoef<-rmvnorm(nsim, coeff, summary(model)$cov.unscaled)
    newdata<-model$family$linkinv(model.matrix(model)%*%t(rcoef))
    # make new datasets from fitted values of gee or gam:
    newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
    for(i in 1:nsim){
      newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata[,i], d = disp) 
    }
  }else{
    # make new datasets from fitted values of gee or gam:
    newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
    for(i in 1:nsim){
      newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata, d = disp) 
    }  
  }
  return(newdata2ndsim)
}

generateDistribData2_bin<-function(nsim, model, newdata, betasamp=FALSE, coeff=NULL, dist.func='quasipoisson', disp=1){
  
  if(dist.func=='quasipoisson'){
    if(betasamp==TRUE){
      rcoef<-rmvnorm(nsim, coeff, summary(model)$cov.unscaled)
      newdata<-model$family$linkinv(model.matrix(model)%*%t(rcoef))
      # make new datasets from fitted values of gee or gam:
      newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
      for(i in 1:nsim){
        newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata[,i], d = disp) 
      }
    }else{
      # make new datasets from fitted values of gee or gam:
      newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
      for(i in 1:nsim){
        newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata, d = disp) 
      }  
    }
  }
  if(dist.func=='binomial'){
    if(betasamp==TRUE){
      rcoef<-rmvnorm(nsim, coeff, summary(model)$cov.unscaled)
      newdata<-model$family$linkinv(model.matrix(model)%*%t(rcoef))
      # make new datasets from fitted values of gee or gam:
      newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
      for(i in 1:nsim){
        newdata2ndsim[,i]<-rbinom(n = nrow(newdata), size = 1, prob = newdata[,i]) 
      }
    }else{
      # make new datasets from fitted values of gee or gam:
      newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
      for(i in 1:nsim){
        newdata2ndsim[,i]<-rbinom(n = nrow(newdata), size=1, prob = newdata) 
      }  
    }
  }
  return(newdata2ndsim)
}
# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

rpois.od<-function(n, lambda, d=1){
  if(d[1]==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

acffunc_dat<-function (block, data) 
{
  blocktab <- table(block)
  acfmat <- matrix(NA, length(unique(block)), max(blocktab))
  for (i in 1:length(unique(block))) {
    print(i)
    corr <- as.vector(acf(data[which(block == unique(block)[i])], plot = F, lag.max = max(blocktab))$acf)
    if(length(which(is.na(corr)))>0){
      corr<-c(1, rep(0.001, (length(corr)-1)))
    }
    acfmat[i, 1:length(corr)] <- corr
  }
  return(list(acfmat = acfmat, blocktab = blocktab))
}


acffunc_dat2<-function (block, data) 
{
  blocktab <- table(block)
  acfmat <- matrix(NA, length(unique(block)), max(blocktab))
  datademean<-data-mean(data)
  for (i in 1:length(unique(block))) {
    print(i)
    corr <- as.vector(acf(datademean[which(block == unique(block)[i])], plot = F, lag.max = max(blocktab), demean = FALSE)$acf)
    #if(length(which(is.na(corr)))>0){
    #  corr<-c(1, rep(0.001, (length(corr)-1)))
    #}
    acfmat[i, 1:length(corr)] <- corr
  }
  return(list(acfmat = acfmat, blocktab = blocktab))
}

# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~
generateIC<-function(data, corrs, block, newdata, nsim){

require(Hmisc)
require(Matrix)
bids<-unique(data[,block])
#nsim=ncol(newdata)
numRep=nsim # number of draws to be taken

for(iter in 1:length(bids)){
  cat(iter , '\n')
  
  tempid<-which(data[,block]==bids[iter])
  
  vars<-t(newdata[tempid,1:nsim])
  
  numVar=length(tempid)  # number of variables to consider
  
  bcorr<-na.omit(corrs$acfmat[iter, ])
  if(bcorr[1]!=1){
    bcorr[1]<-1
  }
  sigma<-matrix(NA, nrow=length(bcorr), ncol=length(bcorr))
  for(i in 1:nrow(sigma)){
    sigma[i,(i:ncol(sigma))]<-bcorr[1:(ncol(sigma)-(i-1))]
  }
  sigma[lower.tri(sigma)]<-t(sigma)[lower.tri(sigma)]
  
  entryR=qnorm((1:numRep)/(numRep+1))
  absDiff=numVar*numVar-numVar
  for (j in 1:200) { #Pick a "good" r matrix - note: 100 is different to numrep
    thisR=matrix(NA, nrow=numRep, ncol=numVar)
    for (i in 1:numVar) {
      thisR[,i]=entryR[sample(1:numRep,numRep)]
      #thisR=cbind(thisR,entryR[sample(1:numRep,numRep)])
    }
    if (j<2) {R=thisR}
    thisAbsDiff=sum(abs(rcorr_sub(thisR)$r))
    if (thisAbsDiff<absDiff) {
      absDiff=thisAbsDiff
      R=thisR
    }
    #cat(thisAbsDiff,', ' , absDiff, '\n')
  }
  
  P=t(chol(sigma)) 
  R_star=R%*%t(P)
  T=rcorr_sub(R)
  # NB T=rcorr(R_star)#from dodgey Haas code
  
  Q=try(t(chol(T$r)), silent=TRUE)
  if(class(Q)=='try-error'){
    Q=t(chol(nearPD(T$r)$mat))  
  }
  
  S=P%*%solve(Q)
  Rb_star=R%*%t(S)
  #NB Rb_star=R_star%*%t(S)#from dodgey Haas code
  repVars=matrix(nrow=numRep,ncol=numVar)
  #RrepVars=matrix(nrow=numRep,ncol=numVar)
  if(numVar>1){
    iin<-apply(vars, 2, sort)
    outs<-apply(Rb_star, 2, order)
  }else{
    iin<-vars
    outs<-t(t(order(Rb_star)))
  }
              
  for (i in 1:numVar) {
    #repVars[order(Rb_star[,i]),i]=sort(vars[,i])
    repVars[outs[,i],i]<-iin[,i]
   # RrepVars[order(R_star[,i]),i]=sort(vars[,i])  
  }
  #rcorr(repVars)$r
  #rcorr(RrepVars)$r
  
  if(iter==1){
    totalrepVars<-NULL
    #totalRrepVars<-NULL
  }
  
  totalrepVars <- rbind(totalrepVars, t(repVars))
  #totalRrepVars <- rbind(totalRrepVars, t(RrepVars))
}
return(totalrepVars)
}


# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~
generateIC_meancor<-function(data, corrs, block, newdata, nsim){
  
  # corrs is the output from acf with max block length
  require(Hmisc)
  require(Matrix)
  bids<-unique(data[,block])
  #nsim=ncol(newdata)
  numRep=nsim # number of draws to be taken
  
  for(iter in 1:length(bids)){
    cat(iter , '\n')
    
    tempid<-which(data[,block]==bids[iter])
    
    vars<-t(newdata[tempid,1:nsim])
    
    numVar=length(tempid)  # number of variables to consider
    
    bcorr<-corrs[1:length(tempid)]
    # if(bcorr[1]!=1){
    #   bcorr[1]<-1
    # }
    sigma<-matrix(NA, nrow=length(bcorr), ncol=length(bcorr))
    if(length(bcorr)==1){
      sigma<-matrix(1, 1, 1)
    }else{
      for(i in 1:nrow(sigma)){
        sigma[i,(i:ncol(sigma))]<-bcorr[1:(ncol(sigma)-(i-1))]
      }
    }
    sigma[lower.tri(sigma)]<-t(sigma)[lower.tri(sigma)]
    
    entryR=qnorm((1:numRep)/(numRep+1))
    absDiff=numVar*numVar-numVar
    for (j in 1:200) { #Pick a "good" r matrix - note: 100 is different to numrep
      thisR=matrix(NA, nrow=numRep, ncol=numVar)
      for (i in 1:numVar) {
        thisR[,i]=entryR[sample(1:numRep,numRep)]
        #thisR=cbind(thisR,entryR[sample(1:numRep,numRep)])
      }
      if (j<2) {R=thisR}
      thisAbsDiff=sum(abs(rcorr_sub(thisR)$r))
      if (thisAbsDiff<absDiff) {
        absDiff=thisAbsDiff
        R=thisR
      }
      #cat(thisAbsDiff,', ' , absDiff, '\n')
    }
    
    P=t(chol(sigma)) 
    R_star=R%*%t(P)
    T=rcorr_sub(R)
    # NB T=rcorr(R_star)#from dodgey Haas code
    
    Q=try(t(chol(T$r)), silent=TRUE)
    if(class(Q)=='try-error'){
      Q=t(chol(nearPD(T$r)$mat))  
    }
    
    S=P%*%solve(Q)
    Rb_star=R%*%t(S)
    #NB Rb_star=R_star%*%t(S)#from dodgey Haas code
    repVars=matrix(nrow=numRep,ncol=numVar)
    #RrepVars=matrix(nrow=numRep,ncol=numVar)
    if(numVar>1){
      iin<-apply(vars, 2, sort)
      outs<-apply(Rb_star, 2, order)
    }else{
      iin<-vars
      outs<-t(t(order(Rb_star)))
    }
    
    for (i in 1:numVar) {
      #repVars[order(Rb_star[,i]),i]=sort(vars[,i])
      repVars[outs[,i],i]<-iin[,i]
      # RrepVars[order(R_star[,i]),i]=sort(vars[,i])  
    }
    #rcorr(repVars)$r
    #rcorr(RrepVars)$r
    
    if(iter==1){
      totalrepVars<-NULL
      #totalRrepVars<-NULL
    }
    
    totalrepVars <- rbind(totalrepVars, t(repVars))
    #totalRrepVars <- rbind(totalRrepVars, t(RrepVars))
  }
  return(totalrepVars)
}

# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~
rcorr_sub<-function (x, y, type = c("pearson", "spearman")) 
{
  type <- match.arg(type)
  if (!missing(y)) 
    x <- cbind(x, y)
  x[is.na(x)] <- 1e+50
  storage.mode(x) <- "double"
  p <- as.integer(ncol(x))
  if (p < 1) 
    stop("must have >1 column")
  n <- as.integer(nrow(x))
  if (n < 5) 
    stop("must have >4 observations")
  h <- .Fortran("rcorr", x, n, p, itype = as.integer(1 + (type == 
                                                            "spearman")), hmatrix = double(p * p), npair = integer(p * 
                                                                                                                     p), double(n), double(n), double(n), double(n), double(n), 
                integer(n), PACKAGE = "Hmisc")
  npair <- matrix(h$npair, ncol = p)
  h <- matrix(h$hmatrix, ncol = p)
  h[h > 1e+49] <- NA
  nam <- dimnames(x)[[2]]
  dimnames(h) <- list(nam, nam)
  dimnames(npair) <- list(nam, nam)
  # P <- matrix(2 * (1 - pt(abs(h) * sqrt(npair - 2)/sqrt(1 - 
  #                                                         h * h), npair - 2)), ncol = p)
  # P[abs(h) == 1] <- 0
  # diag(P) <- NA
  # dimnames(P) <- list(nam, nam)
  structure(list(r = h, n = npair, P = c()), class = "rcorr_sub")
}

# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

getBetaCoverage<-function(betas, ses, df, coefs){
  uppers<- betas+qt(0.025, df=df, lower.tail = FALSE)*ses
  lowers<- betas-qt(0.025, df=df, lower.tail = FALSE)*ses
  inside<-ifelse(coefs>lowers & coefs<uppers,1,0)
return(inside)
}


# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

simfunc_SW<- function(nsim, model, impactcoeff, datatype='Independent', simdata, correlations=NULL){
  
  # find the correction for the se's and also the correction for the variance-covariance matrix (to some extent these are the same)
  data<-model$data
  # find eventphase correction
  fitglm_noev<-glm(model$formula, data=data, family=quasipoisson)
  fitglm<-update(fitglm_noev, .~. + eventphase)
  fitgee<-update(model, .~. + eventphase, data=data)
  correc_imp<-summary(fitgee)$coefficients[,2]/summary(fitglm)$coefficients[,2]
  #correc_imp
  
  
  # generate independent data
  if(is.null(simdata)){
    newdata2ndsim<-generateDistribData(nsim, fitgee, coeff = coef(fitgee))
    if(datatype=='Correlated'){
      # generate correlated data
      cat('Generating Correlated Data...\n')
      newdata2ndsim<-generateIC(data, correlations, "blockid2", newdata2ndsim, nsim=nsim)
    }
  }else{
    newdata2ndsim<-simdata
  }
  
  # make storage outputs
  betas=matrix(NA, length(coef(fitgee)), nsim)
  ses_glm=ses_gee=ses_cimp=matrix(NA, length(coef(fitgee)), nsim)
  disp=vector(length=nsim)
  inside_glm=inside_gee= inside_corr=matrix(NA, nsim, length(coef(fitgee)))
  pvalsimp<-matrix(NA, nsim, 3)
  
  # re-fit model (gam) with new response vectors. 
  for(i in 1:nsim){
    print(i)
    data$simsimresp<-newdata2ndsim[,i]
    
    # glm impact fit
    simfit_glmimp<-update(fitglm, simsimresp ~., data=data)
    pvalsimp[i,1]<-summary(simfit_glmimp)$coefficients[6,4]
    betas[,i]<-summary(simfit_glmimp)$coefficients[,1]
    ses_glm[,i]<-summary(simfit_glmimp)$coefficients[,2]
    inside_glm[i,]<-getBetaCoverage(betas[,i], ses_glm[,i], summary(simfit_glmimp)$df[2], c(coef(model), impactcoeff))
    
    # gee impact fit (should be no different as no correlation actually present)
    simfit_geeimp<-update(fitgee, simsimresp ~., data=data)
    pvalsimp[i,2]<-summary(simfit_geeimp)$coefficients[6,4]
    ses_gee[,i]<-summary(simfit_geeimp)$coefficients[,2]
    disp[i]<-as.numeric(summary(simfit_geeimp)$dispersion[1])
    inside_gee[i,]<-getBetaCoverage(betas[,i], ses_gee[,i], summary(simfit_geeimp)$df[2], c(coef(model), impactcoeff))
    
    # find the corrected se's and see if coverage is better.
    ses_cimp[,i]<-summary(simfit_glmimp)$coefficients[,2]*correc_imp
    inside_corr[i,]<-getBetaCoverage(betas[,i], ses_cimp[,i], summary(simfit_geeimp)$df[2], c(coef(model), impactcoeff))
    
    # take glm output and use corrected standard errors to make new teststat
    teststat<-summary(simfit_glmimp)$coefficients[6,1]/ses_cimp[6,i]
    # and p-value
    pvalsimp[i,3]<-2*pt(-abs(teststat), summary(simfit_glmimp)$df[2])
  }
  
  glm.out<-list(sderr=ses_glm, sderr.corrected=ses_cimp, betas=betas, beta.coverage=inside_glm)
  gee.out<-list(sderr=ses_gee, betas=betas, disp=disp, beta.coverage=inside_glm, beta.coverage.corrected=inside_corr)
  
  output<-list(correc_imp=correc_imp,glm.out=glm.out, gee.out=gee.out, p.evphase=pvalsimp, newdata2ndsim=newdata2ndsim)
  return(output)
}



simfunc_SW2<- function(nsim, model, initialcoeff, datatype='Independent', simdata, correlations=NULL){
  
  # find the correction for the se's and also the correction for the variance-covariance matrix (to some extent these are the same)
  data<-model$data
  # find eventphase correction
  fitglm_noev<-glm(model$formula, data=data, family=quasipoisson)
  fitglm<-update(fitglm_noev, .~. + eventphase)
  fitgee<-update(model, .~. + eventphase, data=data)
  correc_imp<-summary(fitgee)$coefficients[,2]/summary(fitglm)$coefficients[,2]
  #correc_imp
  
  
  # generate independent data
  if(is.null(simdata)){
    newdata2ndsim<-generateDistribData(nsim, fitgee, coeff = coef(fitgee))
    if(datatype=='Correlated'){
      # generate correlated data
      cat('Generating Correlated Data...\n')
      newdata2ndsim<-generateIC(data, correlations, "blockid2", newdata2ndsim, nsim=nsim)
    }
  }else{
    newdata2ndsim<-simdata
  }
  
  # make storage outputs
  betas=matrix(NA, length(coef(fitgee)), nsim)
  ses_glm=ses_gee=ses_cimp=matrix(NA, length(coef(fitgee)), nsim)
  disp=vector(length=nsim)
  inside_glm=inside_gee= inside_corr=matrix(NA, nsim, length(coef(fitgee)))
  pvalsimp<-matrix(NA, nsim, 3)
  
  # re-fit model (gam) with new response vectors. 
  for(i in 1:nsim){
    print(i)
    data$simsimresp<-newdata2ndsim[,i]
    
    # glm impact fit
    simfit_glmimp<-update(fitglm, simsimresp ~., data=data)
    pvalsimp[i,1]<-summary(simfit_glmimp)$coefficients[6,4]
    betas[,i]<-summary(simfit_glmimp)$coefficients[,1]
    ses_glm[,i]<-summary(simfit_glmimp)$coefficients[,2]
    inside_glm[i,]<-getBetaCoverage(betas[,i], ses_glm[,i], summary(simfit_glmimp)$df[2], initialcoeff)
    
    # gee impact fit (should be no different as no correlation actually present)
    simfit_geeimp<-update(fitgee, simsimresp ~., data=data)
    pvalsimp[i,2]<-summary(simfit_geeimp)$coefficients[6,4]
    ses_gee[,i]<-summary(simfit_geeimp)$coefficients[,2]
    disp[i]<-as.numeric(summary(simfit_geeimp)$dispersion[1])
    inside_gee[i,]<-getBetaCoverage(betas[,i], ses_gee[,i], summary(simfit_geeimp)$df[2], initialcoeff)
    
    # find the corrected se's and see if coverage is better.
    ses_cimp[,i]<-summary(simfit_glmimp)$coefficients[,2]*correc_imp
    inside_corr[i,]<-getBetaCoverage(betas[,i], ses_cimp[,i], summary(simfit_geeimp)$df[2], initialcoeff)
    
    # take glm output and use corrected standard errors to make new teststat
    teststat<-summary(simfit_glmimp)$coefficients[6,1]/ses_cimp[6,i]
    # and p-value
    pvalsimp[i,3]<-2*pt(-abs(teststat), summary(simfit_glmimp)$df[2])
  }
  
  glm.out<-list(sderr=ses_glm, sderr.corrected=ses_cimp, betas=betas, beta.coverage=inside_glm)
  gee.out<-list(sderr=ses_gee, betas=betas, disp=disp, beta.coverage=inside_glm, beta.coverage.corrected=inside_corr)
  
  output<-list(correc_imp=correc_imp,glm.out=glm.out, gee.out=gee.out, p.evphase=pvalsimp, newdata2ndsim=newdata2ndsim)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~
generateDistribData_sim<-function(nsim, model, coeff, data, disp, dist.func=NULL){
  rcoef<-rmvnorm(nsim, coeff, summary(model)$cov.unscaled)
  newdata<-model$family$linkinv(model.matrix(model)%*%t(rcoef))
  
  # make new datasets from fitted values of gee or gam:
  newdata2ndsim<-matrix(NA, nrow=nrow(newdata), ncol=nsim)
  
  if(dist.func=='quasipoisson'){
  for(i in 1:nsim){
    newdata2ndsim[,i]<-rpois.od(n = nrow(newdata), lambda = newdata[,i], d = disp) 
  }
  }
  
  if(dist.func=='binomial'){
      for(i in 1:nsim){
        newdata2ndsim[,i]<-rbinom(n = nrow(newdata), size=1, prob = newdata) 
      }  
    }
  return(newdata2ndsim)
}






simfunc_SW_disp<- function(nsim, model, impactcoeff, datatype='Independent', simdata, correlations=NULL, disp){
  
  # find the correction for the se's and also the correction for the variance-covariance matrix (to some extent these are the same)
  data<-model$data
  # find eventphase correction
  fitglm_noev<-glm(model$formula, data=data, family=quasipoisson)
  fitglm<-update(fitglm_noev, .~. + eventphase)
  fitgee<-update(model, .~. + eventphase, data=data)
  correc_imp<-summary(fitgee)$coefficients[,2]/summary(fitglm)$coefficients[,2]
  #correc_imp
  
  
  # generate independent data
  if(is.null(simdata)){
    newdata2ndsim<-generateDistribData_sim(nsim, fitgee, coeff = coef(fitgee), disp = disp, dist.func='quasipoisson')
    if(datatype=='Correlated'){
      # generate correlated data
      cat('Generating Correlated Data...\n')
      newdata2ndsim<-generateIC(data, correlations, "blockid2", newdata2ndsim, nsim=nsim)
    }
  }else{
    newdata2ndsim<-simdata
  }
  
  # make storage outputs
  betas=matrix(NA, length(coef(fitgee)), nsim)
  ses_glm=ses_gee=ses_cimp=matrix(NA, length(coef(fitgee)), nsim)
  disp=vector(length=nsim)
  inside_glm=inside_gee= inside_corr=matrix(NA, nsim, length(coef(fitgee)))
  pvalsimp<-matrix(NA, nsim, 3)
  
  # re-fit model (gam) with new response vectors. 
  for(i in 1:nsim){
    print(i)
    data$simsimresp<-newdata2ndsim[,i]
    
    # glm impact fit
    simfit_glmimp<-update(fitglm, simsimresp ~., data=data)
    pvalsimp[i,1]<-summary(simfit_glmimp)$coefficients[6,4]
    betas[,i]<-summary(simfit_glmimp)$coefficients[,1]
    ses_glm[,i]<-summary(simfit_glmimp)$coefficients[,2]
    inside_glm[i,]<-getBetaCoverage(betas[,i], ses_glm[,i], summary(simfit_glmimp)$df[2], c(coef(model), impactcoeff))
    
    # gee impact fit (should be no different as no correlation actually present)
    simfit_geeimp<-update(fitgee, simsimresp ~., data=data)
    pvalsimp[i,2]<-summary(simfit_geeimp)$coefficients[6,4]
    ses_gee[,i]<-summary(simfit_geeimp)$coefficients[,2]
    disp[i]<-as.numeric(summary(simfit_geeimp)$dispersion[1])
    inside_gee[i,]<-getBetaCoverage(betas[,i], ses_gee[,i], summary(simfit_geeimp)$df[2], c(coef(model), impactcoeff))
    
    # find the corrected se's and see if coverage is better.
    ses_cimp[,i]<-summary(simfit_glmimp)$coefficients[,2]*correc_imp
    inside_corr[i,]<-getBetaCoverage(betas[,i], ses_cimp[,i], summary(simfit_geeimp)$df[2], c(coef(model), impactcoeff))
    
    # take glm output and use corrected standard errors to make new teststat
    teststat<-summary(simfit_glmimp)$coefficients[6,1]/ses_cimp[6,i]
    # and p-value
    pvalsimp[i,3]<-2*pt(-abs(teststat), summary(simfit_glmimp)$df[2])
  }
  
  glm.out<-list(sderr=ses_glm, sderr.corrected=ses_cimp, betas=betas, beta.coverage=inside_glm)
  gee.out<-list(sderr=ses_gee, betas=betas, disp=disp, beta.coverage=inside_glm, beta.coverage.corrected=inside_corr)
  
  output<-list(correc_imp=correc_imp,glm.out=glm.out, gee.out=gee.out, p.evphase=pvalsimp, newdata2ndsim=newdata2ndsim)
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~
simfunc_SW_disp_slim<- function(nsim, model, impactcoeff, correc_imp, datatype='Independent', simdata, correlations=NULL, disp, geemodel=NULL){
  
  data<-model$data
  # generate independent data
  if(is.null(simdata)){
    newdata2ndsim<-generateDistribData_sim(nsim, model, coeff = coef(model), disp = disp, dist.func = 'quasipoisson')
    if(datatype=='Correlated'){
      # generate correlated data
      cat('Generating Correlated Data...\n')
      newdata2ndsim<-generateIC(data, correlations, "blockid2", newdata2ndsim, nsim=nsim)
    }
  }else{
    newdata2ndsim<-simdata
  }
  
  # make storage outputs
  #betas=matrix(NA, length(coef(model)), nsim)
  ses_cimp=matrix(NA, length(coef(model)), nsim)
  #disp=vector(length=nsim)
  #inside_glm=inside_gee= inside_corr=matrix(NA, nsim, length(coef(model)))
  pvalsimp<-matrix(NA, nsim, 2)
  
  
  # re-fit model (gam) with new response vectors. 
  for(i in 1:nsim){
    print(i)
    data$simsimresp<-newdata2ndsim[,i]
    
    impcoefid<-length(coef(model))
    # glm impact fit
    simfit_glmimp<-update(model, simsimresp ~., data=data)
    pvalsimp[i,1]<-summary(simfit_glmimp)$coefficients[impcoefid,4]
    
    # # gee impact fit (should be no different as no correlation actually present)
    if(datatype=='Correlated'){
      simfit_geeimp<-update(geemodel, simsimresp ~., data=data)
      pvalsimp[i,2]<-summary(simfit_geeimp)$coefficients[impcoefid,4]
    }else{
      # find the corrected se's and see if coverage is better.
      ses_cimp[,i]<-summary(simfit_glmimp)$coefficients[,2]*correc_imp
      # take glm output and use corrected standard errors to make new teststat
      teststat<-summary(simfit_glmimp)$coefficients[impcoefid,1]/ses_cimp[impcoefid,i]
      # and p-value
      pvalsimp[i,2]<-2*pt(-abs(teststat), summary(simfit_glmimp)$df[2])  
    }
  }
  
  #glm.out<-list(sderr.corrected=ses_cimp)
  #gee.out<-list(sderr=ses_gee, betas=betas, disp=disp, beta.coverage=inside_glm, beta.coverage.corrected=inside_corr)
  
  output<-list(p.evphase=pvalsimp)
  return(output)
}

# ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~
simfunc_SW_bin_slim<- function(nsim, model, impactcoeff, correc_imp, datatype='Independent', simdata, correlations=NULL, geemodel=NULL){
  
  data<-model$data
  # generate independent data
  if(is.null(simdata)){
    newdata2ndsim<- generateDistribData_sim(nsim, model, coeff=coef(model), dist.func = 'binomial')
    if(datatype=='Correlated'){
      # generate correlated data
      cat('Generating Correlated Data...\n')
      newdata2ndsim<-generateIC(data, correlations, "blockid2", newdata2ndsim, nsim=nsim)
    }
  }else{
    newdata2ndsim<-simdata
  }
  
  # make storage outputs
  #betas=matrix(NA, length(coef(model)), nsim)
  ses_cimp=matrix(NA, length(coef(model)), nsim)
  #disp=vector(length=nsim)
  #inside_glm=inside_gee= inside_corr=matrix(NA, nsim, length(coef(model)))
  pvalsimp<-matrix(NA, nsim, 2)
  
  
  # re-fit model (gam) with new response vectors. 
  for(i in 1:nsim){
    print(i)
    data$simsimresp<-newdata2ndsim[,i]
    
    impcoefid<-length(coef(model))
    # glm impact fit
    simfit_glmimp<-update(model, simsimresp ~., data=data)
    pvalsimp[i,1]<-summary(simfit_glmimp)$coefficients[impcoefid,4]
    
    # # gee impact fit (should be no different as no correlation actually present)
    if(datatype=='Correlated'){
      simfit_geeimp<-update(geemodel, simsimresp ~., data=data)
      pvalsimp[i,2]<-summary(simfit_geeimp)$coefficients[impcoefid,4]
    }else{
      # find the corrected se's and see if coverage is better.
      ses_cimp[,i]<-summary(simfit_glmimp)$coefficients[,2]*correc_imp
      # take glm output and use corrected standard errors to make new teststat
      teststat<-summary(simfit_glmimp)$coefficients[impcoefid,1]/ses_cimp[impcoefid,i]
      # and p-value
      pvalsimp[i,2]<-2*pt(-abs(teststat), summary(simfit_glmimp)$df[2])  
    }
  }
  
  #glm.out<-list(sderr.corrected=ses_cimp)
  #gee.out<-list(sderr=ses_gee, betas=betas, disp=disp, beta.coverage=inside_glm, beta.coverage.corrected=inside_corr)
  
  output<-list(p.evphase=pvalsimp)
  return(output)
}