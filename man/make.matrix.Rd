\name{make.matrix}
\alias{make.matrix}
\title{Utility function used by \code{get.adjacency.matrix()}}
\description{Utility function used by \code{get.adjacency.matrix()}}
\usage{
make.matrix(vertices, interactions);
}
\arguments{
  \item{vertices}{Comma separated list of nodes}
  \item{interactions}{Comma separated list of edges}
}
\value{Returns adjacency matrix}
\author{Syed Haider}
\examples{
x1 <- make.matrix("a,b,c", "a:b,b:c");
}
\keyword{Networks}
