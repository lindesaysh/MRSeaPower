# make boundary

quilt.plot(lincs$x.pos, lincs$y.pos, lincs$depth, asp=1, nrow=80, ncol=80)
bnd<-locator()
bnd<-data.frame(bnd)
write.csv(bnd, file='data/lincsboundary.csv', row.names=F)

# make grid

x<-seq(min(bnd[,1]), max(bnd[,1]), length=100)
y<-seq(min(bnd[,2]), max(bnd[,2]), length=100)

grid<-expand.grid(x, y)
require(splancs)

grid2<-grid[inout(pts = grid, poly = bnd),]
plot(grid2)


