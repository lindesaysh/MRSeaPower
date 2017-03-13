#' Function to make non-regular data into a raster
#'
#' @description This function is useful for when the aim is to calculate the differences between two surfaces which may not have exactly the same grid.
#'
#' @param ncell.y number of raster cells required in the y dimension.
#' @param xyzdata data frame containing c(x.pos, y.pos, truth)
#' @param z.name Character variable denoting the name to be given to the z covariate
#'
#' @author LAS Scott-Hayward, University of St Andrews


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
