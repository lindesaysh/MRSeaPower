#' Function plot the mean, lower 2.5\% and upper 97.5\% of predicted animal counts before and after the event.
#'
#' @param power.object power analysis object of class gamMRSea.power
#' @param returndata logical stating whether to return the plotting data
#' @param cellsize vector of length two denoting the size of the gridcells (width, height).  If unspecified, geom_raster is used for plotting rather than geom_tile.
#'
#' @examples
#'
#' # Load a power object that was created with a 20\% sitewide decline
#' # 50 simulations were run
#' data(nysted.power.oc)
#'
#' plot.preds(nysted.power.oc)
#'
#' @author LAS Scott-Hayward, University of St Andrews
#'
#' @export
#'
plot.preds<-function(power.object, returndata=FALSE, cellsize){

  data<-power.object$bootpreds

  if(length(grep('raster', search()))==1){
    detach(grep('raster', search()), character.only=TRUE)
  }
  library(dplyr)

  a<-select(data, x.pos, y.pos, eventphase, mean)%>%
    mutate(type='mean')
  b<-select(data, x.pos, y.pos, eventphase, X2.5.)%>%
    mutate(type='Lower_2.5')
  b<-rename(b, mean=X2.5.)
  c<-select(data, x.pos, y.pos, eventphase, X97.5.)%>%
    mutate(type='Upper_97.5')
  c<-rename(c, mean=X97.5.)

  plotdata<-data.frame(rbind(a, b, c))

  if(returndata==TRUE){
    return(plotdata)
  }

  data("mypalette")
  breaks<-c(0, 0.2, 0.4, 0.6, 0.8, 1, 10)
  #p<-ggplot(plotdata) + geom_raster(aes( x.pos , y.pos , fill = mean ) )

  if(is.null(cellsize)){
    p<-ggplot(plotdata) + geom_raster(aes(x.pos , y.pos , fill = mean))
  }else{
    p<-ggplot(plotdata) + geom_tile(aes(x.pos , y.pos , fill = mean, width=cellsize[1], height=cellsize[2]))
  }

  p + scale_fill_gradientn(colours=mypalette, values = breaks, space = "Lab", na.value = "grey50", guide = "colourbar", name='Mean Prediction') + theme_bw() + facet_wrap(~type+eventphase, nrow=3) + coord_equal()

}
