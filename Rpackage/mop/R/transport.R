transport <- function(p1, p2, cost, oType=26, lambda=-1){

   n1 <- length(p1)
   n2 <- length(p2)
   
   res <- .Call("transport", as.double(p1), n1, as.double(p2), n2,
       as.double(cost), as.integer(oType), as.double(lambda));

  sol <- list(plan=res[[1]], cost = res[[2]])

}

