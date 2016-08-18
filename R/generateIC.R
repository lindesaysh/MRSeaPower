generateIC<-function(data, corrs, panels, newdata, nsim, dots=TRUE){

  require(Hmisc)
  require(Matrix)
  bids<-unique(data[,panels])
  #nsim=ncol(newdata)
  numRep=nsim # number of draws to be taken

  for(iter in 1:length(bids)){
    if(dots){if((iter/100)%%1 == 0){cat(iter, '\n')}else{cat('.')}}

    if(iter==1){
      totalrepVars<-NULL
      #totalRrepVars<-NULL
    }

    tempid<-which(data[,panels]==bids[iter])

    vars<-t(newdata[tempid,1:nsim])

    numVar=length(tempid)  # number of variables to consider

    if(is.matrix(corrs)){
	bcorr<-na.omit(corrs[iter, ])
	}else{
	bcorr<-na.omit(corrs)
	}

  if(length(bcorr)==0){
    repVars<-vars
  }else{

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
} #end else
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
