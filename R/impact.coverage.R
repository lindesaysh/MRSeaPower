#'
#'  Get the coverage for the impact coefficient
#'
#'  @param truebeta
#'  @param betacis
#'
#'  @example
#'
#'  @export
#'
#'  @author Lindesay Scott-Hayward
#'

impact.coverage<-function(truebeta, betacis){
  length(which(betacis[,1]<truebeta  & betacis[,2]>truebeta))/nrow(betacis)
}
