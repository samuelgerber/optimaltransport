\name{multiscale.transport}
\alias{mop}
\alias{multiscale.transport}
\alias{multiscale.transport.solve}
\alias{multiscale.transport.create.lp}
\alias{multiscale.transport.setup.default}
\alias{multiscale.transport.add.expand.neighborhood.strategy}
\alias{multiscale.transport.add.refine.neighborhood.strategy}
\alias{multiscale.transport.add.potential.neighborhood.strategy}
\alias{multiscale.transport.set.propagation.strategy.1}
\alias{multiscale.transport.set.propagation.strategy.2}
\alias{multiscale.transport.set.max.neighborhood.size}
\alias{multiscale.transport.create.neighborhood.propagation.strategy}
\alias{multiscale.transport.create.max.entropy.propagation.strategy}
\alias{multiscale.transport.create.sinkhorn.propagation.strategy}
\title{Multiscale Transport }
\description{
Compute optimal transport using a multiscale linear programming approach.
}


\usage{


multiscale.transport.solve(gmra1,  gmra2, scale1=-1, scale2=-1 , p=1, 
    w1=1, w2=1, matchScale=TRUE, multiscaleCost=FALSE, multiscaleSolution=FALSE,
    dType=1, nType=1, scaleMass = TRUE)



multiscale.transport.create.lp( oType=31, lambda=-1)
multiscale.transport.setup.default()
multiscale.transport.add.expand.neighborhood.strategy(trp, expandFactor,
    tolerance, iterations=1, maxAdd=10000000)

multiscale.transport.add.refine.neighborhood.strategy(trp, expandFactor,
    tolerance, iterations=1, maxAdd=10000000)

multiscale.transport.add.potential.neighborhood.strategy(trp, threshold,
    tolerance, iterations=1, sort=TRUE, expandPotential=FALSE, maxAdd=10000000)

multiscale.transport.set.propagation.strategy.1(trp, prop)
multiscale.transport.set.propagation.strategy.2(trp, prop)
multiscale.transport.set.max.neighborhood.size(trp, n)

multiscale.transport.create.neighborhood.propagation.strategy(factor = 0)
multiscale.transport.create.randomized.neighborhood.propagation.strategy(nRandom, factor = 0)
multiscale.transport.create.capacity.propagation.strategy(k = 3, factor = 0)
multiscale.transport.create.iterated.capacity.propagation.strategy(nIter = 3, factor = 0)
multiscale.transport.create.max.entropy.propagation.strategy()
multiscale.transport.create.sinkhorn.propagation.strategy(lambda, tolerance,
    threshold, iterations=100)

}

\arguments{
  \item{scale}{Start optimization from scale. 0 = finest scale only (i.e. no
      multiscale optimization), negative all scales, otherwise n scales up from
    finest sacle}
  \item{scale1}{See scale}
  \item{scale2}{See scale}
  \item{stpPct}{Stop neighborhood expansion if improvement of previous solution was
    less than stpPct of previous solution.}
  \item{rFactor}{Multiplicative factor to expand neighborhood when using
    nieghborhood refinment or expand stratgey.}
  \item{p}{p-th Wasssertstein distance}
  \item{nType}{Neighborhood type, 0=absolut, 1 = relative}
  \item{sType}{Neighborhood expansion strategy 0 = Expand, 1 = Refine,
    2=Potential, 3=Fixed, 4=Expand-Potential, 5=Expand-Correct,
    6=Refine-Potential, 7=Refine-Correct. Expand and
      refine use geomtric balls about current source and target locations to add
  new paths to consider in the solution at the next finer scale. Potential uses the
  reduced cost to find all possible paths that need to be added at the next
  scale. Expand-Potential first expands with balls and the uses the potential
  stratgey. Expand-Correct uses the epxand strategy until the last scale at
  which the potential stratgey is employed to find an exact solution.
  Refine-Potential, Refine-Correct are analogous but use the refine strategy
  instead of the expand strategy. }
  \item{dType}{Cost function type. 1 = Euclidean metric, 2 = L1 metric, 3 =
    Squared Euclidean Metric, 4 = Gaussian 2-Wasserstein per Node}
  \item{oType}{Optimizer type: 
    0-2: GLPK primal simplex, 
    11: MOSEK free simplex,
    12: MOSEK network simplex, 
    13: MOSEK free, 
    14: MOSEK interior point, 
    15: MOSEK concurrent,
    21: CPLEX Auto,
    22: CPLEX Barrier,
    23: CPLEX Dual,
    24: CPLEX Primal,
    25: CPLEX Sifting,
    26: CPLEX Network,
    31: Lemon Network Simplex
    Default is lemon
  }
  \item{propFactor}{Neighborhood size (ball = propFactor * node radius) to
                    consider when propating solution from one scale to the next. 
                    0 selects children of current solution only, i.e. arcs with 
                    non-zero flow.}
  \item{nRefinmentIterations}{Number of refinment strategy iterations. 
                              If set to -1 the refinment is iterated unitl convergence.}
  \item{gmra1}{GMRA handle for source distribution (see package gmra)}
  \item{gmra2}{GMRA handle for target distribution (see package gmra)}
  \item{wStein}{Setup gmra such that wassserstein node distance can be used,
                i.e. compute covariance structure in each node}
  \item{nTrials}{Number of randomized runs}
  \item{nExclusions}{Number of exclusion runs}
  \item{factor}{Number of additional paths for refinement in the duality
                strategy to add as a multiplicatve factor of the number of 
                points at each scale}  
  \item{dType}{Cost function type. 1 = Euclidean metric, 2 = L1 metric, 3 =
               Squared Euclidean Metric, 4 = Gaussian 2-Wasserstein per Node}
  \item{lambda}{If lambda is > 0 do a fuzzy match allwoing a mismatch at each
                each target of +/- lambda * mass of node. Only works with CPLEX 
                lp-solvers.}
  \item{k}{Constrain capacity of paths to 1/k times the minimum of the mass of 
           source or target}
  \item{nIter}{Number of capcity constrain iterations}
  \item{scaleMass}{Scale mass of source and taget to sum to one. Otherwise the
                   method solves for an imbalanced mass transport for which some 
                   mass  can remain stationary (the amount that differs between 
                   source and target masses). } }

\value{
A multiscale transport plan for multiscale.transport and
multiscale.transport.ipca.id. It consists of the variables for each scale
(either as a vector or list):\cr
cost - optimal transport cost\cr
nvars - number of paths (variables) in the optimal solution\cr
costsTo - total cost contribution of each target point (sum = optimal transport)\cr
costsFrom - total cost contribution od each source point (sum = optimal tansport cost)\cr
potFrom - dual variable or potential for source points \cr 
potTo - dual variable or potential for source points \cr
map - tranpsort plan, matrix with 3 columns: source index, target index and mass\cr
nTotalPaths - number of total paths considered\cr
from -  source points\cr
to - target points  \cr
fromMass - mass for each source point  \cr 
toMass - mass of each target point\cr
}

\author{
  Samuel Gerber
}
\seealso{
  multiscale.transport.plot
}
\references{

}
\examples{

#create example data  
phi <- runif(1000)*2*pi
X1<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)
X1[,1]=X1[,1]*0.5;
X1[,2]=X1[,2]*2;

phi <- runif(1653)*2*pi
X2<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)
X2[,1]=X2[,1]*3;
X2[,2]=X2[,2]*0.5;

#create multiscale decompsotions
library(gmra)
gmra1 = gmra.create.ipca(X=X1, eps=0, d=2, maxKids=2)
gmra2 = gmra.create.ipca(X=X2, eps=0, d=2, maxKids=2)

#setup and solve multiscale lp
library(mop)
trp.lp <- multiscale.transport.create.lp(oType=30)
icprop <- multiscale.transport.create.iterated.capacity.propagation.strategy(1, 0)
multiscale.transport.set.propagation.strategy.1(trp.lp, icprop);

time1 <- system.time( 
    trp1 <- multiscale.transport.solve(trp.lp, gmra1, gmra2, p = 2, nType=0, dType=1) )

#add nieghborhood expansion
multiscale.transport.add.expand.neighborhood.strategy(trp.lp, 1 ) 
time2 <- system.time( 
    trp2 <- multiscale.transport.solve(trp.lp, gmra1, gmra2, p = 2, nType=0, dType=1) )

#solve optimal
time3 <- system.time( 
    trp3 <- multiscale.transport.solve(trp.lp, gmra1, gmra2, p = 2, nType=0, dType=1, scale1=0, scale2=0) )

trp1$cost
trp2$cost
trp3$cost
time1
time2
time3

}
\keyword{multiscale,gmra,neighbor}
