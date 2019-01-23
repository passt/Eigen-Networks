hexpress <- function(xy, exprmat, feats, pcomp, colix) {
  require(grid)
  require(hexbin)
  require(marray)
  # modified hexagon bin function from hexbin package
  source('./3rdparty/gridhexagons2.R')
  
  # xy-boundaries
  bndLo = c(-22,  -9, -23, -11)
  bndUp = c( 22,  12,   9,  11) # NG
  
  # minimum data per bin
  mincnt <- 5
  
  # remove points outside some boundary.
  is.na(xy[,1]) <- xy[,1] < bndLo[pcomp[1]] | xy[,1] > bndUp[pcomp[1]]
  is.na(xy[,2]) <- xy[,2] < bndLo[pcomp[2]] | xy[,2] > bndUp[pcomp[2]]
  
  
  if (feats == 'density') { 
    # Calculate density
    hb12 <- hexbin::hexbin(xy, IDs=TRUE, xbins=40,
                           xbnds=c(bndLo[pcomp[1]], bndUp[pcomp[1]]),
                           ybnds=c(bndLo[pcomp[2]], bndUp[pcomp[2]]))
    
    par(mar=c(2,2,.5,.5), bty='n')
    P <- plot(hb12, type='n', legend=F, xaxt='n', yaxt='n', clip='on', main='Density')
    P$plot.vp@hexVp.off$xscale=c(bndLo[pcomp[1]], bndUp[pcomp[1]])
    P$plot.vp@hexVp.off$yscale=c(bndLo[pcomp[2]], bndUp[pcomp[2]])
    pushViewport(P$plot.vp@hexVp.off)
    grid.hexagons2(hb12, use.count=T, style= "colorscale",
                   colramp=function(n){maPalette(low='#EEEEEE', high='#000000', k = n)}, mincnt=mincnt)
    da <- xy[hb12@cID %in% hb12@cell[hb12@count < mincnt],]
    grid.xaxis(); grid.yaxis()
    grid.points(da[,1], da[,2], pch='.')
    popViewport()
  }
  else if (feats == 'overall'){
    # Calculate total expression level
    a <- rowSums(exprmat[,colix])
    
    # Visualise as hexbin
    hb12 <- hexbin::hexbin(xy, IDs=TRUE)
    hb12@cAtt <- hexTapply(hbin = hb12, FUN = median, dat = a)
    par(mar=c(2,2,.5,.5))
    P <- plot(hb12, type='n', legend=F, xlab=paste0("PC", pcomp[1]), ylab=paste0("PC", pcomp[2]), main='Overall expression')
    pushHexport(P$plot.vp)
    grid.hexagons2(hb12, use.count=F, cell.at=NULL,  style= "colorscale",
                   colramp=function(n){maPalette(low='blue', high='green', k = n)})
    popViewport()    
  }
  else {
    nodeix   <- which(colnames(exprmat) == feats)
    
    hb12 <- hexbin::hexbin(xy, IDs=TRUE)
    hb12@cAtt <- hexTapply(hbin = hb12, FUN = mean, dat = exprmat[,nodeix])
    par(mar=c(2,2,.5,.5))  
    P <- plot(hb12, type='n', legend=F, xlab=paste0("PC", pcomp[1]), ylab=paste0("PC", pcomp[2]), main=feats)
    pushHexport(P$plot.vp)
    grid.hexagons2(hb12, use.count=F, cell.at=NULL,  style= "colorscale", mincnt=0, colramp=function(n){maPalette(low='blue', high='green', k = n)})
    popViewport()
  }
  
}