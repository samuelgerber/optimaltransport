\name{multiscale.sinkorn.transport}
\alias{sinkorn.transport}
\alias{multiscale.sinkhorn.transport}
\title{ Multiscale Sinkhorn Transport }
\description{
Compute optimal transport with entropy penalty (Sinkhorn distances) either
  direct or a sparse approach using a multiscale method.
}


\usage{
sinkhorn.transport(mu, nu, D, lambda=10, maxIter=1000, tol=1e-5)

multiscale.sinkhorn.transport.ipca(gmra1, gmra2, scale1=-1, scale2=-1, p=1, w1
      = 1, w2 = 1, dType=1, nType=0, lambda=10, tolerance = 1e-5,
      maxIter=10000, threshold = 0, maxPathsPerNode = 10,
      multiscale.cost=FALSE)

multiscale.sinkhorn.transport.plot.map(mst, index, plotMap = TRUE, colX1 =
    t( c(0,0,0) ), colX2 = t( c(1,0,0) ), colMap = t( c(0,0,0) ), cex=1, xlab=expression(x[1]),
    ylab=expression(x[2]), mapAlpha = 1, pointAlpha = 1, X1 = mst$from[[index]],
    X2 = mst$to[[index]], asp=1, add=FALSE)
}

\arguments{
  \item{mu}{Source distribution}
  \item{nu}{Target distirbution}
  \item{D}{Cost matrix, i.e. C[i, j] cost of transporting mass from mu[i] to
    nu[0] }
  \item{gmra1}{GMRA handle (see gmra package)}
  \item{gmra2}{GMRA handle (see gmra package)}
  \item{scale1}{Number of scale to include of GMRA approxmiation of X1. 0 =
    finest scale only, negative all scales, otherwise n scales up from finest
      sacle}
  \item{scale2}{Number of scale to include of  GMRA approxmiation of X2. 0 =
    finest scale only, negative all scales, otherwise n scales up from finest
      scale}
  \item{p}{p-th Wasssertstein distance}
  \item{nType}{Neighborhood type, 0=absolut, 1 = relative}
  \item{dType}{Cost function type. 1 = Euclidean metric, 2 = L1 metric, 3 =
    Squared Euclidean Metric, 4 = Gaussian Wasserstein}
  \item{lambda}{Entropy penalty weight, the smaller the larger}
  \item{maxIter}{ Maximum number of iterations}
  \item{tol}{Minimal relative change in transport cost for stopping}
  \item{maxPathsPerNode}{When moving from one scale to the next allow at most
    this many paths per node to be propagted. The transport plan that is
      propgated is pruned to the maxPathsPerNode most mass bearing paths per
      node.}
  \item{threshold}{Paths with less than threshold times maximal mass posible for
  this path are not included in the propagation to the next scale.}
  \item{mst}{multiscale transport plan strcuture}
  \item{index}{scale of the multiscale transport plan to plot}
  \item{plotMap}{draw lines for the transport plan assignments}
  \item{colX1}{a color for the source data points}
  \item{colX2}{a color for the targert data points}
  \item{colMap}{color for the transport plan}
  \item{cex}{magnification factor for drawn points}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{mapAlpha}{alpha factor for drawing transport plan}
  \item{pointsAlpha}{alpha factor for drawing points}
  \item{wStein}{Setup gmra such that wassserstein node distance can be used,
    i.e. compute covariance structure in each node}
}
\value{
A multiscale transport plan structure for multiscale.sinkhorn.transport.ipca and
  the sinkhorn distance and runtime for sinkorn.transport.
}

\author{
  Samuel Gerber
}
\seealso{
}
\references{
Sinkhorn Distances: Lightspeed Computation of Optimal Transportation
Cuturi, Marco
Advances in Neural Information Processing Systems, 2013
}
\examples{

phi <- runif(1000)*2*pi
X1<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

phi <- runif(5003)*2*pi
X2<- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)
X2[,1]=X2[,1]*3;
X2[,2]=X2[,2]*0.5;

library(gmra)
gmra1 = gmra.create.ipca(X=X1, eps=0, d=2, maxKids=2)
gmra2 = gmra.create.ipca(X=X2, eps=0, d=2, maxKids=2)

trp <- multiscale.sinkhorn.transport(gmra1, gmra2, p = 2 );

plot(1:length(trp$cost), trp$cost, xlab="scale", ylab="transport cost / sinkhorn
    distance" ) 

multiscale.sinkhorn.transport.plot.map(trp, length(trp$cost), mapAlpha=0.025) 


}
\keyword{multiscale,gmra,neighbor}
