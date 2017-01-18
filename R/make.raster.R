#' Function to make non-regular data into a raster to calculate the differences.
#'
#' @param ncell.y
#' @param before
#' @param after
#' @param legend.name
#'


make.raster<-function(ncell.y, xyzdata, z.name = 'response'){
  require(raster)

  spdf1<- SpatialPointsDataFrame( data.frame( x = xyzdata$x.pos , y = xyzdata$y.pos ) ,
                                  data = data.frame( z = xyzdata$truth) )
  e <- extent( spdf1 )
  ratio <- ( e@xmax - e@xmin ) / ( e@ymax - e@ymin )
  r <- raster( nrows = ncell.y , ncols = floor( ncell.y * ratio ) , ext = extent(spdf1) )
  rf <- rasterize( spdf1 , r , field = "z" , fun= mean)
  rdf1 <- data.frame( rasterToPoints( rf ))

  names(rdf1)[3]<-z.name

  return(rdf1)
}
