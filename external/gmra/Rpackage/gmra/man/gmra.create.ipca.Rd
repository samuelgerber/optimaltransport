\name{gmra.create.ipca}
\alias{gmra.create.ipca}
\title{Create Iterated PCA Tree}
\description{
Create hierarchical decomposition by recursive splitting, eitehr aong principal
  components or axis aligned.
}


\usage{
gmra.create.ipca( X, eps, dim=3, maxKids=3, threshold = -1, split = 2, splitDir = 1,
    stop=4, nRuns=1, minPoints=1 )
}
\arguments{
  \item{X}{Data set. Each row is a data point}
  \item{eps}{Accuracy of GMRA approximation as defined in stop}
  \item{stop}{Stop if: 1 = node MSE / total MSE  < eps, 2 = node
  R^2 < eps, 3 = node MSE < eps, 4: node radius < eps, 5 = node Radius / total
    Radius < eps, 6: node radius * node mass < eps}
  \item{split}{1 = split at mean, 2 = split at mid point of projected points, 3 =
  split at mean + random varition, 4 = split at mid point + random variation}
  \item{splitDir}{1 = split across principal components, 2 = split aacross
    randomized principal components, 3 = Axis aligned }
  \item{dim}{Maximal dimensionality of each GMRA node}
  \item{maxKids}{Allow at most 2^maxKids children per node}
  \item{threshold}{Threshold for number of principal components retained. If < 0 each
    node has d pc's otherwise if t < 1 it retains as many as required to capyure
      t percent of the variation but at most d }
  \item{nRuns}{Create multiple trees. Useful for creating a set of randomized
      trees.}
  \item{minPoints}{Minimum number of points in a node. Additional stoppinf criterium}
}

\examples{

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

library(gmra)
gmra = gmra.create.ipca(arc, eps=0, dim=2, maxKids=2, stop=4)


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
