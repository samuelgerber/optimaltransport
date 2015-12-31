













multiscale.transport.plot.map <- function(mst, index, plotMap = TRUE, colX1 =
    t( c(0,0,0) ), colX2 = t( c(1,0,0) ), colMap = t( c(0,0,0) ), cex=1, xlab=expression(x[1]),
    ylab=expression(x[2]), mapAlpha = 1, pointAlpha = 1, X1 = mst$from[[index]],
    X2 = mst$to[[index]], asp=1, colorByPotential = FALSE, add=FALSE, lwd=2,
    cex.axis=1, cex.lab=1, arrows=FALSE, arrow.length=0.1, arrow.angle=15){

  library(RColorBrewer)
  if(index > length(mst$cost) ){
    index = length(mst$cost) 
  }

  if(plotMap){
    map = mst$map[[index]]
  }
  else{
    map=NULL
  }

  if(colorByPotential){
    pX1 = mst$potFrom[[index]]
    pX2 = mst$potTo[[index]]

    mi = min(c(pX1, pX2))
    pX1 = pX1 - mi
    pX2 = pX2 - mi
    ma = max(c(pX1, pX2))
    pX1 = pX1 / ma
    pX2 = pX2 / ma
    ramp = colorRamp( brewer.pal(n=11, name="PuOr") ) 
    colX1 = ramp(pX1)/255 
    colX2 = ramp(pX2)/255 
  }

  if(pointAlpha < 0){
    alpha = min(1, -pointAlpha)
    pointAlpha1 = alpha
    pointAlpha2 = alpha
  }
  else{
    z = max(c(mst$fromMass[[index]], mst$toMass[[index]]) )
    pointAlpha1 = mst$fromMass[[index]] / z * pointAlpha
    pointAlpha2 =  mst$toMass[[index]] / z * pointAlpha
    pointAlpha1[pointAlpha1>1] = 1
    pointAlpha2[pointAlpha2>1] = 1
       
  }

  if(mapAlpha < 0){
    mapAlpha = -mapAlpha
    mapAlpha = min(1, mapAlpha)  
  }
  else{
    mapAlpha = mst$map[[index]][,3] / max(mst$map[[index]][,3]) * mapAlpha
    mapAlpha[mapAlpha>1] = 1  
  }

  xlim = range(c(X1[,1], X2[,1]))
  ylim = range(c(X1[,2], X2[,2]))

  if( !is.null(map) ){ 
    if(!add){
        plot( NULL, pch=19, cex=cex, col = rgb(colX1, alpha=pointAlpha1), xlim=xlim,
      ylim=ylim, xlab=xlab, ylab=ylab, bty="n", asp=asp, cex.lab=cex.lab,
      cex.axis=cex.axis )
    }
    if(arrows){
     arrows( x0 = X1[map[,1], 1], y0 = X1[map[,1], 2], x1 = X2[map[,2], 1], y1 =
         X2[map[,2], 2], col=rgb(colMap, alpha=mapAlpha ), lwd=lwd,
         angle=arrow.angle,
         code=2, length=arrow.length )
    }
    else{
      segments( x0 = X1[map[,1], 1], y0 = X1[map[,1], 2], x1 = X2[map[,2], 1], y1 =
          X2[map[,2], 2], col=rgb(colMap, alpha=mapAlpha ), lwd=lwd )
    }
    add=TRUE;

  }

  if(add){
    points( X1, pch=19, cex=cex, col = rgb(colX1, alpha=pointAlpha1) )
  }
  else{
    plot( X1, pch=19, cex=cex, col = rgb(colX1, alpha=pointAlpha1), xlim=xlim,
      ylim=ylim, xlab=xlab, ylab=ylab, bty="n", asp=asp, cex.lab=cex.lab,
      cex.axis=cex.axis )
  }
  points( X2, pch=19, cex=cex, col = rgb(colX2, alpha=pointAlpha2) )
  
}













multiscale.transport.plot.multiscale.map <- function(mst, index, plotMap = TRUE, colX1 =
    t( c(0,0,0) ), colX2 = t( c(1,0,0) ), colMap = t( c(0,0,0) ), cex=1, xlab=expression(x[1]),
    ylab=expression(x[2]), mapAlpha = 1, pointAlpha = 1, asp=1, colorByPotential = FALSE, add=FALSE){

  library(RColorBrewer)  
    
  if(index > length(mst$cost) ){
    index = length(mst$cost) 
  }

  if(plotMap){
    map = mst$map[[index]]
  }
  else{
    map=NULL
  }

  X1 = mst$from[[index]]
  X2 = mst$to[[index]]

  if(colorByPotential){
    pX1 = mst$potFrom[[index]]
    pX2 = mst$potTo[[index]]

    mi = min(c(pX1, pX2))
    pX1 = pX1 - mi
    pX2 = pX2 - mi
    ma = max(c(pX1, pX2))
    pX1 = pX1 / ma
    pX2 = pX2 / ma
    ramp = colorRamp( brewer.pal(n=11, name="PuOr") ) 
    colX1 = ramp(pX1)/255 
    colX2 = ramp(pX2)/255 
  }

  if(pointAlpha < 0){
    alpha = min(1, -pointAlpha)
    pointAlpha1 = alpha
    pointAlpha2 = alpha
  }
  else{
    z = max(c(mst$fromMass[[index]], mst$toMass[[index]]) )
    pointAlpha1 = mst$fromMass[[index]] / z * pointAlpha
    pointAlpha2 =  mst$toMass[[index]] / z * pointAlpha
    pointAlpha1[pointAlpha1>1] = 1
    pointAlpha2[pointAlpha2>1] = 1
       
  }

  if(mapAlpha < 0){
    mapAlpha = -mapAlpha
    mapAlpha = min(1, mapAlpha)  
  }
  else{
    mapAlpha = mst$map[[index]][,3] / max(mst$map[[index]][,3]) * mapAlpha
    mapAlpha[mapAlpha>1] = 1  
  }

  xlim = range(c(X1[,1], X2[,1]))
  ylim = range(c(X1[,2], X2[,2]))
  if(add){
    points( X1, pch=19, cex=cex, col = rgb(colX1, alpha=pointAlpha1) )
  }
  else{
    plot( X1, pch=19, cex=cex, col = rgb(colX1, alpha=pointAlpha1), xlim=xlim,
      ylim=ylim, xlab=xlab, ylab=ylab, bty="n", asp=asp )
  }
  points( X2, pch=19, cex=cex, col = rgb(colX2, alpha=pointAlpha2) )
  
  if( !is.null(map) ){
    XS = mst$from[[index]][mst$mFrom[index, ], ]
    XE = mst$to[[index]][mst$mTo[index, ], ]
    
    D = list()
    Dsum = matrix(0, nrow=ncol(mst$mFrom), ncol=ncol(XS) )
    for(i in 1:index){ 
      D[[i]] = mst$to[[i]][mst$mTo[i,], ] - mst$from[[i]][mst$mFrom[i,], ] - Dsum
      Dsum = Dsum + D[[i]]
    }
    
    
    for(i in index:1){
      D[[i]] = D[[i]]/2
      segments( x0 = XS[, 1], y0 = XS[, 2], x1 = (XS+D[[i]])[, 1], y1 =
          (XS+D[[i]])[, 2], col=rgb(colMap, alpha=mapAlpha ) )
      segments( x0 = XE[, 1], y0 = XE[, 2], x1 = (XE-D[[i]])[, 1], y1 =
          (XE-D[[i]])[, 2], col=rgb(colMap, alpha=mapAlpha ) )
      XS = XS+D[[i]]
      XE = XE-D[[i]]
    }
  }

}










multiscale.transport.plot.binnned.map <- function(trp, bins, xlab="x",
    ylab="y", col="black", arrow.length=0.1, arrow.angle=15, alpha=1, lwd=2,
    add = FALSE, cex.axis=1, cex.lab=1, cex=1, useTransparancy=FALSE,
    useCost=FALSE){

  Q = multiscale.transport.bin(trp, bins)
 
  dx = Q$dx[[1]]
  dy = Q$dx[[2]]

  xloc = Q$loc[[1]]
  yloc = Q$loc[[2]]

  if(useCost){
    ws=Q$cost/max(Q$cost)
  }
  else{
    ws = Q$w/max(Q$w)
  }


  if(!add){
    plot(NA, xlim=range(bins[[1]]), ylim=range(bins[[2]]), xlab=xlab, ylab=ylab,
        cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, bty="n")
  }
  for(i in 1:Q$dims[1]){
    for(j in 1:Q$dims[2]){
      if(useTransparancy){
        arrows( x0 = xloc[i], y0 = yloc[j], x1 = xloc[i]+dx[i, j], y1 = yloc[j]
          +dy[i, j], col=rgb( t(col2rgb(col)/255), alpha=0.05+0.95*ws[i,j]), lwd=lwd, angle=arrow.angle,
          code=2, length=arrow.length )
      }
      else{
        arrows( x0 = xloc[i], y0 = yloc[j], x1 = xloc[i]+dx[i, j], y1 = yloc[j]
          +dy[i, j], col=rgb( t(col2rgb(col)/255), alpha=0.7), lwd=1+ws[i,j]*lwd, angle=arrow.angle,
          code=2, length=arrow.length )
      }
    }
  }

}









multiscale.transport.plot.binnned.map.3d <- function(trp, bins, xlab="x",
    ylab="y", zlab="z", col.from="black", col.to="red", arrow.length=0.1, arrow.angle=15, alpha=1, lwd=2,
    add = FALSE, cex.axis=1, cex.lab=1, cex=1, useTransparancy=FALSE,
    useCost=FALSE){

  library(rgl)
  
    
  Q = multiscale.transport.bin(trp, bins)
 
  dx = Q$dx[[1]]
  dy = Q$dx[[2]]
  dz = Q$dx[[3]]

  xloc = Q$loc[[1]]
  yloc = Q$loc[[2]]
  zloc = Q$loc[[3]]

  if(useCost){
    ws=Q$cost/max(Q$cost)
  }
  else{
    ws = Q$w/max(Q$w)
  }



  if(!add){
    xlim=range(bins[[1]])
    ylim=range(bins[[2]])
    zlim=range(bins[[3]])
    plot3d(cbind(xlim,ylim,zlim), xlim=xlim, ylim=ylim, zlim=zlim, xlab=xlab, ylab=ylab,
        zlab =zlab, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, bty="n",
        col="#00000000")
  }


  for(i in 1:Q$dims[1]){
    for(j in 1:Q$dims[2]){
      for(k in 1:Q$dims[3]){
        x = c(xloc[i], yloc[j], zloc[k])
        V = rbind(x, x+ c(dx[i, j, k], dy[i, j, k],dz[i, j, k]) )


        if(useTransparancy){
          color = c( rgb( t(col2rgb(col.from)/255), alpha=0.05+0.95*ws[i,j,k]),
                   rgb( t(col2rgb(col.to)/255), alpha=0.05+0.95*ws[i,j,k]) )
          lwidth = lwd;
        }
        else{
          color = c( rgb( t(col2rgb(col.from)/255), alpha=0.7),
                   rgb( t(col2rgb(col.to)/255), alpha=0.7) )
          lwidth = lwd=1+ws[i,j,k]*lwd 
        }
        segments3d( V, col=color, lwd=lwidth)
      }

    }
  }

}



multiscale.transport.bin <- function(trp, bins){

 
  l = length(trp$cost)
  plan = trp$map[[l]]
  from = trp$from[[l]]
  to = trp$to[[l]]

  
  n <- list()
  loc <- list()
  cuts <- list()
  dims <- c()
  Ind = matrix(nrow=length(bins), ncol=nrow(plan))
  for(i in 1:length(bins)){
    n[[i]] = length(bins[[i]])
    loc[[i]] = bins[[i]][1:(n[[i]]-1)] + (bins[[i]][2:n[[i]]] -
      bins[[i]][1:(n[[i]]-1)]) / 2 
    Ind[i, ] = cut(from[ plan[, 1], i] , bins[[i]])
    dims <- c(dims, length(bins[[i]])-1)

  }


  w = array(0, dims)
  dx = list()
  delta = to[ plan[, 2], ] - from[ plan[,1], ]
  for(i in 1:ncol(delta)){
    dx[[i]] = w;
  }

  for( i in 1:nrow(plan) ){
    ind = t(Ind[, i])
    w[ind] = w[ind]  + plan[i, 3]
    for(j in 1:length(dx)){
      dx[[j]][ind] = dx[[j]][ind] + plan[i, 3]*delta[i, j]
    }
  }


  for( j in 1:length(dx) ){
    dx[[j]] = dx[[j]]/w
  }


  list(w=w, dx=dx, loc=loc, n=n, dims=dims) 

}


