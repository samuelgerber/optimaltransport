




multiscale.parse.result <- function(res, multiscaleSolution){
   lTo = list()
   lFrom = list()
   pTo = list()
   pFrom = list()
   to = list()
   from = list()
   map = list()
   fromMass = list()
   toMass = list()

   fromIndex = list()
   toIndex = list()
   fromSize = list()
   toSize = list()
   fromRadius = list()
   toRadius = list()

   nTotalPaths = rep( 0, length(res[[1]]) )
   timesRefine = rep( 0, length(res[[1]]) )
   timesLP = rep(0, length(res[[1]]) )

   nVars = 17

   for(i in 1:length(res[[1]]) ){
     ii = i-1
     lFrom[[i]] = res[[ 2 + nVars*ii+1 ]]
     lTo[[i]] = res[[ 2 + nVars*ii+2 ]]
     pFrom[[i]] = res[[ 2 + nVars*ii+3 ]]
     pTo[[i]] = res[[ 2 + nVars*ii+4 ]]
     map[[i]] = res[[ 2 + nVars*ii+5 ]]
     colnames(map[[i]]) = c("from", "to", "mass", "cost")

     nTotalPaths[i] = res[[ 2 + nVars*ii+6 ]]
     from[[i]] = t( res[[ 2 + nVars*ii+7 ]] )
     to[[i]] = t( res[[ 2 + nVars*ii+8 ]] )
     fromMass[[i]] = t( res[[ 2 + nVars*ii+9 ]] )
     toMass[[i]] = t( res[[ 2 + nVars*ii+10 ]] )
     timesLP[i] =  res[[ 2 + nVars*ii+11]][1] 
     timesRefine[i] =  res[[ 2 + nVars*ii+11]][2]

     fromIndex[[i]] = res[[ 2 + nVars*ii+12 ]]
     fromSize[[i]] = res[[ 2 + nVars*ii+13 ]]
     fromRadius[[i]] = res[[ 2 + nVars*ii+14 ]]
     toIndex[[i]] = res[[ 2 + nVars*ii+15 ]]
     toSize[[i]] = res[[ 2 + nVars*ii+16 ]]
     toRadius[[i]] = res[[ 2 + nVars*ii+17 ]]

   }   
   
   mFrom = NULL 
   mTo = NULL 
   mCost = NULL 
   if(multiscaleSolution){
     index = 2+nVars*length(res[[1]])+1
     mFrom = res[[index]] +1
     mTo = res[[index+1]]+1 
     mCost = res[[index+2]] 
   }

   sol = structure( list( cost = res[[1]], nvars = res[[2]], costsTo = lTo, costsFrom =
     lFrom , potFrom = pFrom, potTo = pTo, map = map, nTotalPaths = nTotalPaths,
     from= from, to = to, fromMass = fromMass, toMass = toMass, timesRefine =
     timesRefine, timesLP = timesLP, mCost = mCost, mFrom=mFrom, mTo=mTo,
     fromIndex=fromIndex, toIndex=toIndex, fromSize=fromSize, toSize=toSize,
     fromRadius = fromRadius, toRadius = toRadius),
       class="multiscale.transport" )

  sol

}




multiscale.transport.solve <- function(trp, gmra1, gmra2, scale1=-1, scale2=-1,
    stpPct = 0, p=1, w1=1, w2=1, matchScale=FALSE, multiscaleCost=FALSE,
    multiscaleSolution=FALSE, dType=1, nType=1){

    res <- .Call("multiscaleTransportSolve", trp, gmra1, gmra2, as.integer(scale1),
        as.integer(scale2), as.double(w1), as.integer(length(w1)),
        as.double(w2), as.integer(length(w2)), as.double(p),
        as.integer(matchScale), as.integer(multiscaleCost),
        as.integer(multiscaleSolution), as.integer(dType), as.integer(nType)  )

   res <- multiscale.parse.result(res, multiscaleSolution)
   res$p <- p

}






### setup routines

multiscale.transport.setup.default <- function(lambda=0, oType=31, nIter=1){
   trp <- multiscale.transport.create.lp(oType=oType, lambda=lambda)
   prop <- multiscale.transport.create.iterated.capacity.propagation.strategy(nIter, 0);
   multiscale.transport.set.propagation.strategy.1(trp, prop);
   multiscale.transport.add.expand.neighborhood.strategy(trp, 1);

   trp
}


multiscale.transport.create.lp <- function(oType=31, lambda=-1){
  res <- .Call("createTransportLP", as.integer(oType), as.double(lambda) );
  res
}


#expansion strategies after propagation
multiscale.transport.add.expand.neighborhood.strategy <- function(trp,
    expandFactor, tolerance=0, iterations=1, maxAdd=10000000){
  res <- .Call("addExpandNeighborhoodStrategy", trp, as.double(expandFactor),
      as.double(tolerance), as.integer(iterations), as.integer(maxAdd) )
  res
}


multiscale.transport.add.refine.neighborhood.strategy <- function(trp,
    expandFactor, tolerance=0, iterations=1, maxAdd=10000000){
  res <- .Call("addRefineNeighborhoodStrategy", trp, as.double(expandFactor),
      as.double(tolerance), as.integer(iterations), as.integer(maxAdd) )
  res
}


multiscale.transport.add.potential.neighborhood.strategy <- function(trp,
    threshold, tolerance=0, iterations=1, sort=TRUE, expandPotential=FALSE, maxAdd=10000000){
  res <- .Call("addPotentialNeighborhoodStrategy", trp, as.double(threshold),
      as.double(tolerance), as.integer(iterations), as.integer(sort),
      as.integer(expandPotential), as.integer(maxAdd) )
  res
}


multiscale.transport.set.propagation.strategy.1 <- function(trp, prop){
  .Call("setPropagationStrategy1", trp, prop)
}

multiscale.transport.set.propagation.strategy.2 <- function(trp, prop){
  .Call("setPropagationStrategy2", trp, prop)
}

multiscale.transport.set.max.neighborhood.size <- function(trp, n){
  .Call("setMaxNeighborhoodSize", trp, as.integer(n))
}



#Propagation strategies
multiscale.transport.create.neighborhood.propagation.strategy <- function(factor
    = 0){
  .Call("createNeighborhoodPropagationStrategy", as.double(factor))
}


multiscale.transport.create.capacity.propagation.strategy <- function(k = 3, factor
    = 0){
  .Call("createCapacityPropagationStrategy", as.integer(k), as.double(factor))
}

multiscale.transport.create.iterated.capacity.propagation.strategy <- function(nIter = 3, factor
    = 0){
  .Call("createIteratedCapacityPropagationStrategy", as.integer(nIter), as.double(factor))
}

multiscale.transport.create.randomized.neighborhood.propagation.strategy <-
function(nRandom, factor = 0){
  .Call("createRandomizedNeighborhoodPropagationStrategy", as.integer(nRandom), as.double(factor))
}

multiscale.transport.create.max.entropy.propagation.strategy <- function(){
  .Call("createMaxEntropyPropagationStrategy")
}

multiscale.transport.create.sinkhorn.propagation.strategy <- function(lambda,
    tolerance, threshold, iterations=100){
  .Call("createSinkhornPropagationStrategy", as.double(lambda),
      as.double(tolerance), as.double(threshold), as.integer(iterations) )
}


