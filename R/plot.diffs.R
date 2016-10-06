plot.diffs<-function(power.object, returndata=FALSE){

data<-power.object$bootdifferences

if(length(grep('raster', search()))==1){
  detach(grep('raster', search()), character.only=TRUE)
}

library(dplyr)

a<-select(data, x.pos, y.pos, mean)%>%
  mutate(type='mean')
b<-select(data, x.pos, y.pos, X2.5.)%>%
  mutate(type='Lower_2.5')
b<-rename(b, mean=X2.5.)
c<-select(data, x.pos, y.pos, X97.5.)%>%
  mutate(type='Upper_97.5')
c<-rename(c, mean=X97.5.)

plotdata<-data.frame(rbind(a, b, c))

if(returndata==TRUE){
  return(plotdata)
}

require(RColorBrewer)
breaks<-ticks<-c(0, 0.2, 0.4, 0.6, 0.8, 1, 10, 25)
mypalette<-rev(brewer.pal(length (ticks)-1,"Spectral"))

ggplot(plotdata) + geom_raster(aes( x.pos , y.pos , fill = mean ) ) + scale_fill_gradientn(colours=mypalette, values = breaks, space = "Lab", na.value = "grey50", guide = "colourbar", name='Mean Diff') + theme_bw() + facet_wrap(~type, nrow=1) + coord_equal()

}


plot.preds<-function(power.object, returndata=FALSE){

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

  require(RColorBrewer)
  breaks<-ticks<-c(0, 0.2, 0.4, 0.6, 0.8, 1, 10, 25)
  mypalette<-rev(brewer.pal(length (ticks)-1,"Spectral"))

  ggplot(plotdata) + geom_raster(aes( x.pos , y.pos , fill = mean ) ) + scale_fill_gradientn(colours=mypalette, values = breaks, space = "Lab", na.value = "grey50", guide = "colourbar", name='Mean Prediction') + theme_bw() + facet_wrap(~type+eventphase, nrow=3) + coord_equal()

}
