\name{get.program.defaults}
\alias{get.program.defaults}
\title{A utility function to return the inst/ directory of the installed package and other default settings}
\description{A utility function to return the inst/ directory of the installed package to get the test datasets and other program related data contents}
\usage{
get.program.defaults(networks.database = "default");
}
\arguments{
  \item{networks.database}{Name of the pathway networks database. Default to NCI PID/Reactome/Biocarta i-e "default"}
}
\value{Returns a list of paths to the input directories/files where the contents of this package are installed}
\author{Syed Haider}
\examples{
program.data <- get.program.defaults();
}
\keyword{IO}
