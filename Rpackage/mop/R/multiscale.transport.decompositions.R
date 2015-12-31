






multiscale.transport.decompose <-function(trp, w1, w2){

  n = length(trp$cost)
  from <-c()
  to <-c()

  for(i in 1:n ){
    idx = trp$fromIndex[[i]]
    idy = trp$toIndex[[i]]

    startFrom = c(0, cumsum(trp$fromSize[[i]]) )
    startTo = c(0, cumsum(trp$toSize[[i]]) )
    
    x = rep(0, length(idx) )
    for( j in 1:length(trp$fromSize[[i]]) ){
      fi = j
      ind = which(trp$map[[i]][, 1] == fi)
      d = sum(trp$map[[i]][ind,4] * trp$map[[i]][ind,3]) /   sum(trp$map[[i]][ind,3] )
      x[ idx[ (startFrom[fi]+1):startFrom[fi+1] ] ] = d;
    }
    from = cbind(from, x)

    y = rep(0, length(idy) )
    for( j in 1:length( trp$toSize[[i]] ) ){
      ti = j
      ind = which(trp$map[[i]][, 2] == ti)
      d = sum(trp$map[[i]][ind,4] * trp$map[[i]][ind,3]) /   sum(trp$map[[i]][ind,3] )
      y[ idy[ (startTo[ti]+1):startTo[ti+1] ] ] = d;
    }
    to = cbind(to, y)

  }

  list(from=from, to =to)

}




multiscale.transport.decompose2 <-function(trp, w1, w2){
  library(Matrix)

  n = length(trp$cost)
  D <- list()
  W <- list()
  R <- list()

  for(i in 1:n ){
    idx = trp$fromIndex[[i]]
    idy = trp$toIndex[[i]]

    startFrom = c(0, cumsum(trp$fromSize[[i]]) )
    startTo = c(0, cumsum(trp$toSize[[i]]) )
    
    D1 = matrix(0, nrow=length(idx), ncol=length(idy)) 
    W1 = matrix(0, nrow=length(idx), ncol=length(idy)) 
    R1 = matrix(0, nrow=length(idx), ncol=length(idy)) 
    for( j in 1:nrow(trp$map[[i]]) ){
      fi = trp$map[[i]][j, 1]
      ti = trp$map[[i]][j, 2]
      ix = idx[ (startFrom[fi]+1):startFrom[fi+1] ];
      iy = idy[ (startTo[ti]+1):startTo[ti+1] ]
      w = outer(w1[ix], w2[iy])
      W1[ix, iy] = w/sum(w) * trp$map[[i]][j, 3]
      R1[ix, iy] = trp$fromRadius[[i]][fi] +  trp$toRadius[[i]][ti]                          
      D1[ix, iy] = trp$map[[i]][j, 4]
    }

    D[[i]] = D1
    W[[i]] = W1
    R[[i]] = R1
  }

  list(D=D, W=W, R=R)
}




multiscale.transport.vector.decomposition <-function(trp, X1, X2){
  library(Matrix)

  n = length(trp$cost)
  D <- list()
  W <- list()
  R <- list()


  for(i in n:1 ){
    
    idx = trp$fromIndex[[i]]
    idy = trp$toIndex[[i]]

    startFrom = c(0, cumsum(trp$fromSize[[i]]) )
    startTo = c(0, cumsum(trp$toSize[[i]]) )
   
    delta <- matrix(0, nrow=nrow(X1), ncol=ncol(X1))
    w <- rep(0, nrow(X1))
    r <- rep(0, nrow(X1))
    for( j in 1:nrow(trp$map[[i]]) ){
      fi = trp$map[[i]][j, 1]
      ti = trp$map[[i]][j, 2]
      ix = idx[ (startFrom[fi]+1):startFrom[fi+1] ];
      iy = idy[ (startTo[ti]+1):startTo[ti+1] ]
      wtmp = trp$map[[i]][j, 3]
      
      if(i==n){
        w[ix] = w[ix] + wtmp/(length(ix)*length(iy))
        delta[ix, ] = delta[ix, ]    +  wtmp * (trp$to[[i]][ti, ] - trp$from[[i]][fi, ])                
        r[ix] = r[ix] + wtmp/(length(ix)*length(iy)) * ( trp$fromRadius[[i]][fi] +
                                                               trp$toRadius[[i]][ti] )  
      }
      else{
        wtmp2 = sum(W[[n]][ix])
        w[ix] = w[ix] + wtmp2
        delta[ix, ] = delta[ix, ]    +  sum(wtmp * D[[n]][ix]) 
        r[ix] = r[ix] + wtmp2 * ( trp$fromRadius[[i]][fi] + trp$toRadius[[i]][ti] )  
      }
    }

    D[[i]] = delta/w 
    R[[i]] = r/w
    W[[i]] = w/sum(w) 
  }

  list(D=D, W=W, R=R)
}


