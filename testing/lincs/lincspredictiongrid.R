# make boundary

quilt.plot(lincs$x.pos, lincs$y.pos, lincs$depth, asp=1, nrow=80, ncol=80)
bnd<-locator()
bnd<-data.frame(bnd)
write.csv(bnd, file='data/lincsboundary.csv', row.names=F)

# make grid

bnd<-read.csv('data/lincsboundary.csv')

x<-seq(min(bnd[,1]), max(bnd[,1]), length=100)
y<-seq(min(bnd[,2]), max(bnd[,2]), length=100)

grid<-expand.grid(x, y)
require(splancs)

grid2<-grid[inout(pts = grid, poly = bnd),]
plot(grid2, asp=1)
grid2<-na.omit(grid2)


require(akima)
newdepth<-interpp(x=lincs$x.pos, y=lincs$y.pos, z=lincs$depth, xo=grid2$Var1, yo=grid2$Var2, linear=FALSE, extrap=TRUE)
predgrid<-data.frame(newdepth)
predgrid<-predgrid[predgrid$z>0,]

names(predgrid)<-c('x.pos', 'y.pos', 'depth')

par(mfrow=c(1,2))
quilt.plot(lincs$x.pos, lincs$y.pos, lincs$depth, asp=1, nrow=80, ncol=80)
#points(grid2, pch=20)
quilt.plot(predgrid$x.pos, predgrid$y.pos, predgrid$depth, asp=1, nrow=98, ncol=98)
#head(grid2)


write.csv(predgrid, file='data/lincspreddata.csv', row.names=F)
