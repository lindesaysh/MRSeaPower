
#' Summary function for objects of class gamMRSea.power
#'
#' @param power.object  An object of class gamMRSea.power
#' @param
#'
#'
summary.gamMRSea.power<-function(power.object, null.object=NULL, truebeta){

  nsim.m<-length(power.object$imppvals)
  model.power<-(length(which(power.object$imppvals<=0.05))/nsim.m)*100
  model.coverage<-impact.coverage(truebeta, power.object$betacis)*100

  if(is.null(power.object$Mean.Proportion)){
    proptable<-power.object$Abundance
    modeltype='count'
  }else{
    proptable<-power.object$Mean.Proportion
    modeltype='proportion'
  }

  if(is.null(null.object)){
    null.power=NULL
    null.coverage=NULL
  }else{
    nsim.n<-length(null.object$imppvals)
    null.power<-(length(which(null.object$imppvals<=0.05))/nsim.n)*100
    null.coverage<-impact.coverage(0, null.object$betacis)*100
  }

  ans<-list(power=list(model=model.power, null=null.power), coverage=list(model.coverage=model.coverage, null.coverage=null.coverage), summary=proptable, modeltype=modeltype, truebeta=truebeta)

  class(ans) <- "summary.gamMRSea.power"

  return(ans)
}



print.summary.gamMRSea.power<-function (x, digits = max(3L, getOption("digits") - 3L), ...)
{

  cat('\n ++++ Summary of Power Analysis Output ++++\n')

  if(is.null(x$power$null)){

    cat("\nPower to select 'change' term:\n")
    cat('\n    Under Change (truth = ', format(truebeta, digits) , ') = ', x$power$model, '%\n', sep='')
    cat('    Under no change = Null distribution not specified\n')

    cat("\nCoverage for 'change' coefficient:\n")
    cat('\n    Under change = ', x$coverage$model.coverage, '%\n', sep='')
    cat('    Under no change = Null distribution not specified\n', sep='')

  }else{

    cat("\nPower to select 'change' term:\n")
    cat('\n    Under Change (truth = ', format(truebeta, digits) , ') = ', x$power$model, '%\n', sep='')
    cat('    Under no change (truth = ', 0 ,') = ', x$power$null, '%\n', sep='')

    cat("\nCoverage for 'change' coefficient:\n")
    cat('\n    Under model = ', x$coverage$model.coverage, '%\n', sep='')
    cat('    Under no change = ', x$coverage$null.coverage, '%\n', sep='')
  }

  if(x$modeltype=='count'){
    cat('\nOverall Abundance Summary with 95% Confidence Intervals:\n\n')
    print(x$summary, digits=digits)
  }
  if(x$modeltype=='proportion'){
    cat('\nOverall Mean Proportion Summary with 95% Confidence Intervals:\n\n')
    print(x$summary)
  }

  cat('\nNote: These calculations assume the correct area has been given for each prediction grid cell.\n')

}

