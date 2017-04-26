neighbor.min.flow <- function(X1, X2, weights, lambda= 0){

   m <- nrow(X1)
   n <- ncol(X2)
   
   res <- .Call("neighborMinFlow", as.double(X1), as.double(X2), as.double(weights), as.integer(m),
                             as.integer(n), as.double(lambda) );

  sol <- list(plan=res[[1]], id2s = res[[2]], cost = res[[3]])

}




neighbor.min.flow.plot <- function( nmf, colMap = t( c(0,0,0) ), cex=1,
                                    xlab=expression(x[1]),
                                    ylab=expression(x[2]), mapAlpha = 1, asp=1,
                                    add=FALSE, lwd=2, cex.axis=1, cex.lab=1,
                                    arrows=FALSE, arrow.length=0.1,
                                    arrow.angle=15){ 
  
  id2s = nmf$id2s
  
  xlim = range( id2s[,1] )
  ylim = range( id2s[,2] )
  if( !add ){
    plot( NULL, pch=19, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
          bty="n", asp=asp, cex.lab=cex.lab, cex.axis=cex.axis )
  }



  map = nmf$plan
  
  if(mapAlpha < 0){
    mapAlpha = -mapAlpha
    mapAlpha = min(1, mapAlpha)  
  }
  else{
    mapAlpha = map[,3] / max( map[,3]) * mapAlpha
    mapAlpha[ mapAlpha>1 ] = 1  
  }

  if(arrows){
    arrows( x0 = id2s[map[,1], 1], y0 = id2s[map[,1], 2], x1 = id2s[map[,2], 1], y1 =
            id2s[map[,2], 2], col=rgb(colMap, alpha=mapAlpha ), lwd=lwd,
            angle=arrow.angle, code=2, length=arrow.length )
  }
  else{
    segments( x0 = id2s[map[,1], 1], y0 = id2s[map[,1], 2], x1 = id2s[map[,2], 1], y1 =
              id2s[map[,2], 2], col=rgb(colMap, alpha=mapAlpha ), lwd=lwd )
  }

  
}













neighbor.min.flow3d <- function(X1, X2, lambda= 0){

  dims <- dim(X1)
   
  res <- .Call("neighborMinFlow3d", as.double(X1), as.double(X2), as.integer(dims[1]),
                             as.integer(dims[2]), as.integer(dims[3]), as.double(lambda) );

  sol <- list(plan=res[[1]], id2s = res[[2]], cost = res[[3]])

}






neighbor.min.flow3d.plot <- function( nmf, slice, sliceAxis=3, colMap = t( c(0,0,0) ), cex=1,
                                    xlab=expression(x[1]),
                                    ylab=expression(x[2]), mapAlpha = 1, asp=1,
                                    add=FALSE, lwd=2, cex.axis=1, cex.lab=1,
                                    arrows=FALSE, arrow.length=0.1,
                                    arrow.angle=15){ 
  

  id2s = nmf$id2s
  

  xlim = range( id2s[ ,1] )
  ylim = range( id2s[ ,2] )
  if( !add ){
    plot( NULL, pch=19, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
          bty="n", asp=asp, cex.lab=cex.lab, cex.axis=cex.axis )
  }


  map = nmf$plan
 
   
  index = which( id2s[map[,1], sliceAxis] == slice )
  axisD = setdiff(1:3, sliceAxis)
  a1 = axisD[1]
  a2 = axisD[2]

  if(mapAlpha < 0){
    mapAlpha = -mapAlpha
    mapAlpha = min(1, mapAlpha)  
  }
  else{
    mapAlpha = map[index,3] / max( map[index,3]) * mapAlpha
    mapAlpha[ mapAlpha>1 ] = 1  
  }

  if(arrows){
    arrows( x0 = id2s[map[index, 1], a1], y0 = id2s[map[index,1], a2], 
            x1 = id2s[map[index, 2], a1], y1 = id2s[map[index,2], a2], 
            col=rgb(colMap, alpha=mapAlpha ), lwd=lwd,
            angle=arrow.angle, code=2, length=arrow.length )
  }
  else{
    segments( x0 = id2s[map[index,1], a1], y0 = id2s[map[index,1], a2], 
              x1 = id2s[map[index,2], a1], y1 = id2s[map[index,2], a2], 
              col=rgb(colMap, alpha=mapAlpha ), lwd=lwd )
  }
  
}






neighbor.min.flow3d.plot3d <- function( nmf, colMap = t( c(0,0,0) ), cex=1,
                                    xlab=expression(x[1]),
                                    ylab=expression(x[2]), mapAlpha = 1, asp=1,
                                    add=FALSE, lwd=2, cex.axis=1, cex.lab=1,
                                    arrows=FALSE, arrow.length=0.1,
                                    arrow.angle=15){ 
  

  library(rgl)

  id2s = nmf$id2s
  

  xlim = range( id2s[ ,1] )
  ylim = range( id2s[ ,2] )
  zlim = range( id2s[ ,2] )
  if( !add ){
    plot3d( NULL, xlim=xlim, ylim=ylim, zlim=zlim )
  }


  map = nmf$plan
  
  if(mapAlpha < 0){
    mapAlpha = -mapAlpha
    mapAlpha = min(1, mapAlpha)  
  }
  else{
    mapAlpha = map[,3] / max( map[,3]) * mapAlpha
    mapAlpha[ mapAlpha>1 ] = 1  
  }

  x = c( rbind( id2s[map[,1], 1], id2s[map[,2], 1] ) )
  y = c( rbind( id2s[map[,1], 2], id2s[map[,2], 2] ) )
  z = c( rbind( id2s[map[,1], 3], id2s[map[,2], 3] ) )

  if(arrows){
    arrows3d( x=x, y=y, z=z, col=rgb(colMap), alpha=mapAlpha, lwd=lwd,
            angle=arrow.angle, code=2, length=arrow.length )
  }
  else{
    segments3d( x=x, y=y, z=z, col=rgb(colMap), alpha=mapAlpha, lwd=lwd )
  }
  
}



neighbor.min.flow.interpolate <- function(X1, mf, lambda){
  X1 = X1/sum(X1)
  for(i in 1:nrow(mf$plan) ){
    x = mf$id2s[ mf$plan[i, 1], 1] + 1
    y = mf$id2s[ mf$plan[i, 1], 2] + 1
    xx = mf$id2s[ mf$plan[i, 2], 1] + 1
    yy = mf$id2s[ mf$plan[i, 2], 2] + 1
    w = lambda*mf$plan[i, 3] 
    X1[x,y] = X1[x,y] - w
    X1[xx,yy] = X1[xx,yy] + w
  }
  X1
}
