\name{gmra}
\alias{gmra}
\alias{gmra.root}
\alias{gmra.nodes}
\alias{gmra.node.children}
\alias{gmra.node.parent}
\alias{gmra.node.radii}
\alias{gmra.node.npoints}
\alias{gmra.node.centers}
\alias{gmra.node.points}
\alias{gmra.centers}
\alias{gmra.partition}
\alias{gmra.radii}
\alias{gmra.npoints}
\alias{gmra.neighborhood}
\alias{gmra.density.estimate}
\alias{gmra.prune.min.points.at.scale}
\title{GMRA based data analysis}
\description{
Methods to manipulate geometric multresolution analysis trees.
}


\usage{
  gmra.radii(gmra, scale, dType=1)
  gmra.npoints(gmra, scale)
  gmra.neighborhood(x, gmra, scale, eps, dType=1) 
  gmra.centers(gmra, scale)
  gmra.partition(gmra, scale)
  gmra.root(gmra)
  gmra.nodes(gmra, scale)
  gmra.node.parent(node)
  gmra.node.children(node)
  gmra.node.radii(nodes, dType=1)
  gmra.node.npoints(nodes)
  gmra.node.centers(nodes)
  gmra.node.points(nodes)
  gmra.density.estimate(X, gmras, minPoints, maxPoints, sigma, dType=1)
  gmra.prune.min.points.at.scale(gmra, scale, npoints)
}
\arguments{
  \item{gmras}{List of gmra handles}
  \item{gmra}{gmra handle}
  \item{nodes}{List of node handles}
  \item{node}{node handle}
  \item{scale}{Query at scale}
  \item{eps}{Size of nieghborhood}
  \item{dType}{Distance type to use to compare ndoes, 1: L2, 2: L1, 3: Squared
    L2, 4: Wasserstein (requires that the GMRA tree is decorated by GWT, i.e.
        estimate of covariance) }
  \item{minPoints}{use only nodes with at least minPoints points for estimating
  density}
  \item{maxPoints}{use only nodes with at most maxPoints points for estimating
  density}
  \item{sigma}{Kernel bandwidth for density estimation}
  \item{x}{Query point}
  \item{X}{estimate density at points in X}
  \item{npoints}{Minimal number of points at children of scale to avoid pruning}
}

\examples{

##Noisy circle example
phi <- runif(100)*2*pi
arc <- cbind(cos(phi), sin(phi)) * (1+rnorm(length(phi)) * 0.1)

library(gmra)
gmra = gmra.create.ipca(arc, eps=0, dim=2, maxKids=2, stop=4)


}

\value{
}

\author{
  Samuel Gerber
}
\seealso{
  gmra.create.ipca, gmra.create.ikm
}
\references{

}
\keyword{multiscale,gmra,tree}
