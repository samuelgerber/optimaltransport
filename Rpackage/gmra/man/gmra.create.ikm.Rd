\name{gmra.create.ikm}
\alias{gmra.create.ikm}
\title{Create Iterated K-means Tree}
\description{
Create hierarchical decomposition by recursive k-means.
}


\usage{
gmra.create.ikm( X, eps, nKids=64, threshold = 0.001, maxIter=100, split=2,
    stop=4, nRuns=1, similarity=1, minPoints=1 )
}
\arguments{
  \item{X}{Data set. Each row is a data point}
  \item{eps}{Accuracy of GMRA approximation as defined in stop}
  \item{stop}{Stop if 1: = node MSE / total MSE  < eps, 2: MSE < eps, 3: node
    radius < eps, 4: = node Radius / total Radius < eps}
  \item{nKids}{Maximal or fixed number of kids}
  \item{split}{1 = split into nKids children, 2 = split into at most nKids
    children or radius of children < radius/2 of current, 3 = split into at most
      nKids or radius of children < root Radius / 2^scale of children}
  \item{threshold}{Kmeans stop if realtive improvement of mse is less than threshold}
  \item{maxIter}{Maximum Kmeans iterations}
  \item{nRuns}{Create multiple trees. Useful for creating a set of randomized
      trees.}
  \item{similarity}{1: Run K-means on L2 metric, 2: Assumes points are on unit
    sphere and runs K-means on 1-correlation}
  \item{minPoints}{Minimum number of points in a node. Additional stoppinf criterium}
}

\examples{

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

library(gmra)
gmra = gmra.create.ikm(arc, eps=0)

}

\value{
 a list of GMRA tree handles
}

\author{
  Samuel Gerber
}
\seealso{
  gmra, gmra.create.ikm
}
\references{

}
\keyword{multiscale,gmra,tree}
