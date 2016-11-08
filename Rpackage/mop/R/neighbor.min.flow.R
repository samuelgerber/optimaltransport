neighbor.min.flow <- function(X1, X2){

   m <- nrow(X1)
   n <- ncol(X2)
   
   res <- .Call("neighborMinFlow", as.double(X1), as.double(X2), as.integer(m),
                             as.integer(n) );

  sol <- list(plan=res[[1]], id2s = res[[2]], cost = res[[3]])

}

