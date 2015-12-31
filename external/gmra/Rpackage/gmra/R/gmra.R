
.onLoad <- function(libname, pkgname){
    .Call("gmra_init");
}



gmra.create.ipca <- function(X, eps, dim=3, maxKids=3, threshold = -1, split = 2, splitDir = 1,
    stop=4, nRuns=1, minPoints=1 ){
   n <- ncol(X)
   m <- nrow(X)
   res <- .Call( "gmra_create_ipca", as.double(t(X)), m, n, as.integer(dim),
       as.double(threshold), as.double(eps), as.integer(maxKids), 
       as.integer(split), as.integer(splitDir), as.integer(stop),
       as.integer(nRuns), as.integer(minPoints)  )  
  
   if(nRuns > 1){ 
     ans = list()
       for(i in 1:length(res)){
         ans[[i]] <- structure( res[[i]], class = "gmra.tree", gmra.type="ipca" )
       }
   }
   else{
     ans = structure( res[[1]], class = "gmra.tree", gmra.type="ipca" )
   }

   ans
}


gmra.create.ikm <- function(X, eps, nKids=64, threshold = 0.001, maxIter = 100, 
    split = 2, stop=4, nRuns=1, similarity = 1, nSpatial =0, nCor= 0, minPoints=1){
   n <- ncol(X)
   m <- nrow(X)
   res <- .Call( "gmra_create_ikm", as.double(t(X)), m, n, as.double(threshold),
       as.double(eps), as.integer(nKids), as.integer(maxIter),
       as.integer(split), as.integer(stop), as.integer(nRuns),
       as.integer(similarity), as.integer(nSpatial), as.integer(nCor), as.integer(minPoints) )  
  
   if(nRuns > 1){ 
     ans = list()
       for(i in 1:length(res)){
         ans[[i]] <- structure( res[[i]], class = "gmra.tree", gmra.type="ikm" )
       }
   }
   else{
     ans = structure( res[[1]], class = "gmra.tree", gmra.type="ikm" )
   }

   ans
}







gmra.density.estimate <- function(X, gmras, minPoints, maxPoints, sigma, dType=1){
  
  res <- .Call("gmra_density_estimate", as.double(t(X)), ncol(X), nrow(X),
      as.list(gmras), as.integer(minPoints), as.integer(maxPoints),
      as.double(sigma), as.integer(dType) )
  res 

}


gmra.mean.signal.estimate <- function(signal, X, gmras, minPoints){
  
  res <- .Call("gmra_mean_signal_estimate", as.double(signal), length(signal),
      as.double(t(X)), ncol(X), nrow(X), as.list(gmras), as.integer(minPoints) )
  res 

}


gmra.root <- function(gmra){
  res <- .Call("gmra_get_root",  gmra)
  res  
}

gmra.node.children <- function(node){
  res <- .Call("gmra_node_children",  node)
  res  
}

gmra.node.parent <- function(node){
  res <- .Call("gmra_node_parent",  node)
  res  
}

gmra.nodes <-  function(gmra, scale){
  res <- .Call("gmra_nodes",  gmra, as.integer(scale) )  
  res
}

gmra.partition <-  function(gmra, scale){
  res <- .Call("gmra_partition",  gmra, as.integer(scale) )  
  for(i in 1:length(res) ){
    res[[i]] = res[[i]] +1;
  }
  res
}



gmra.node.points <-  function(nodes){
  res <- .Call("gmra_node_points",  nodes)  
  for(i in 1:length(res) ){
    res[[i]] = res[[i]] +1;
  }
  res
}



gmra.neighborhood <-  function(x, gmra, scale, eps, dType=1){
  res <- .Call("gmra_neighborhood", as.double(x), length(x), gmra,
        as.double(eps), as.integer(scale), as.integer(dType) )  
  t(res)
}


gmra.centers <-  function(gmra, scale){
  res <- .Call("gmra_centers", gmra, as.integer(scale) )  
  t(res)
}

gmra.node.centers <-  function(nodes){
  res <- .Call("gmra_node_centers", nodes )  
  t(res)
}



gmra.npoints <-  function(gmra, scale){
  res <- .Call("gmra_npoints", gmra, as.integer(scale) )  
  res
}

gmra.node.npoints <-  function(nodes){
  res <- .Call("gmra_node_npoints", nodes)  
  res
}

gmra.radii <-  function(gmra, scale, dType=1){
  res <- .Call("gmra_radii", gmra, as.integer(scale), as.integer(dType) )  
  res
}

gmra.node.radii <-  function(nodes, dType=1){
  res <- .Call("gmra_node_radii", nodes, as.integer(dType) )  
  res
}


gmra.prune.min.points.at.scale <- function(gmra, scale, npoints){
  .Call("gmra_prune_min_points_scale", gmra, as.integer(scale),
      as.integer(npoints) )  
}


