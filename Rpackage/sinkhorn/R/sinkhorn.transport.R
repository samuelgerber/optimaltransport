sinkhorn.transport <- function(mu, nu, D, lambda=10, maxIter=1000, tol=1e-5){

  
   if( is.null( nrow(nu) ) ){
     m2 = 1
     n2 = length(nu)
   }else{
     n2 <- nrow(nu)
     m2 <- ncol(nu)
   }

   res <- .Call("sinkhorn_transport", as.double(mu), length(mu),
       as.double(t(nu)), n2, m2, as.double(D), as.double(lambda),
       as.double(tol), as.integer(maxIter) );

   sol = structure( list( cost = res[[1]], time = res[[2]]), class="sinkorn.transport" )

   sol
}



multiscale.sinkhorn.transport <- function(gmra1, gmra2, scale1=-1,
    scale2=-1, p=1, w1 = 1, w2 = 1, dType=1, nType=0, lambda=10, tolerance =
    1e-5, maxIter=10000, threshold = 0, maxPathsPerNode = 10,
    multiscale.cost=FALSE){

   n1 <- ncol(X1)
   m1 <- nrow(X1)
   
   n2 <- ncol(X2)
   m2 <- nrow(X2)
   
   res <- .Call("multiscale_sinkhorn_transport", gmra1, gmra2,
       as.integer(scale1), as.integer(scale2), as.double(w1),
       as.integer(length(w1)), as.double(w2), as.integer(length(w2)),
       as.double(p), as.integer(dType), as.integer(nType), as.double(lambda),
       as.double(tolerance), as.integer(maxIter), as.double(threshold),
       as.integer(maxPathsPerNode), as.integer(multiscale.cost) )

   lTo = list()
   lFrom = list()
   pTo = list()
   pFrom = list()
   to = list()
   from = list()
   map = list()
   fromMass = list()
   toMass = list()

   nTotalPaths = rep(0, length(res[[1]]) )

   for(i in 1:length(res[[1]]) ){
     lFrom[[i]] = res[[ 2 + 10*i-9 ]]
     lTo[[i]] = res[[ 2 + 10*i-8 ]]
     pFrom[[i]] = res[[ 2 + 10*i-7 ]]
     pTo[[i]] = res[[ 2 + 10*i-6 ]]
     map[[i]] = res[[ 2 + 10*i-5 ]]
     nTotalPaths[i] = res[[ 2 + 10*i-4 ]]
     from[[i]] = t( res[[ 2 + 10*i-3 ]] )
     to[[i]] = t( res[[ 2 + 10*i-2 ]] )
     fromMass[[i]] = t( res[[ 2 + 10*i-1 ]] )
     toMass[[i]] = t( res[[ 2 + 10*i ]] )
   }
   sol = structure( list( cost = res[[1]], nvars = res[[2]], costsTo = lTo, costsFrom =
     lFrom , potFrom = pFrom, potTo = pTo, map = map, nTotalPaths = nTotalPaths,
     from= from, to = to, fromMass = fromMass, toMass = toMass),
       class="multiscale.sinkhorn.transport" )

  sol

}


multiscale.sinkhorn.transport.plot.map <- function(mst, index, plotMap = TRUE, colX1 =
    t( c(0,0,0) ), colX2 = t( c(1,0,0) ), colMap = t( c(0,0,0) ), cex=1, xlab=expression(x[1]),
    ylab=expression(x[2]), mapAlpha = 1, pointAlpha = 1, X1 = mst$from[[index]],
    X2 = mst$to[[index]], asp=1, add=FALSE){

  library(RColorBrewer)
  if(plotMap){
    map = mst$map[[index]]
  }
  else{
    map=NULL
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
    segments( x0 = X1[map[,1], 1], y0 = X1[map[,1], 2], x1 = X2[map[,2], 1], y1 =
      X2[map[,2], 2], col=rgb(colMap, alpha=mapAlpha ) )
  }

}



multiscale.sinkhorn.transport.interpolate <- function(mst, index, t){
  from =  mst$from[[index]]
  to =  mst$to[[index]]
  map = mst$map[[index]]

  X = (1-t) * from[map[,1], ] + t * to[map[,2], ]

  list(X=X, w=map[,3])

}




sample.from.array <- function(A, unit.dim = F){
  dims = dim(A)
  d = length(dims)
  
  X <- arrayInd( 1:prod(dims), dims )

  w <- A[X]
  w <- w / sum(w)

  if(unit.dim){
    for(i in 1:d){
      X[,i] = X[,i]/dims[i]
    }
  }
  list( X = X[w>0, ], w = w[w>0] )

}




sample.to.array <-function(X, w, dim = ceiling( apply(X, 2, max)) ){
  A <- array(0, dim=dim)
  
#grid = expand.grid(split(t(matrix(c(1/3, 0, -1/3), ncol=3, nrow=ncol(X))), 1:ncol(X))) 
  for(i in 1:nrow(X)){
    if(w[i] == 0){
      next
    }
#   for(k in 1:nrow(grid)){
#     x =  t( as.integer( floor(X[i, ] +  grid[k, ]) ) ) 
#     A[x] = A[x] + w[i] 
#   }
     x =  t( as.integer( round(X[i, ] ) ) ) 
     A[x] = A[x] + w[i] 
  }
#A = A/nrow(grid)
  A
}

