


multiscale.transport.interpolate <- function(mst, index, t){
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


sample.from.array2 <- function(A, unit.dim = F, n=10000){
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
  X = X[w>0, ]
  w = w[w>0]
  w = w/sum(w)
  m = rmultinom(n=1, prob=w, size=n)
  P = matrix(0, nrow=n, ncol=d)
  index = 1;
  for(i in 1:length(m)){
    if(m[i] > 0){
      for(k in 1:m[i]){
        P[index, ] = X[i, ] + runif(d)-0.5
        index = index+1
      }
    }
  }

 P 
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
