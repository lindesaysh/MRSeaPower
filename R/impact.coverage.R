#'
#'  Get the coverage for the impact coefficient
#'
#'  @param truebeta true coefficient value
#'  @param betacis dataframe of upper and lower intervals, where the number of rows is the number of bootstrap simulations
#'
#'
#'  @export
#'
#'  @author Lindesay Scott-Hayward, university of St Andrews
#'

impact.coverage<-function(truebeta, betacis){
  length(which(betacis[,1]<truebeta  & betacis[,2]>truebeta))/nrow(betacis)
}
