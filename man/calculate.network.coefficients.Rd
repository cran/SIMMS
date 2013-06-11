\name{calculate.network.coefficients}
\alias{calculate.network.coefficients}
\title{Calculate Cox statistics for input dataset}
\description{Function to compute hazard ratios for the genes in pathway-derived networks, by aggregating input datasets into one training cohort. The hazard ratios are computed for each pair by calculating the HR of each gene independently and as an interaction (i.e. y = HR(A) + HR(B) + HR(A:B)}
\usage{
calculate.network.coefficients(
  data.directory = ".", output.directory = ".",
  training.datasets = NULL, data.types = c("mRNA"),
  subnets.file.flattened = NULL, subset = NULL
  );
}
\arguments{
  \item{data.directory}{Path to the directory containing datasets as specified by \code{training.datasets}}
  \item{output.directory}{Path to the output folder where intermediate and results files will be saved}
  \item{training.datasets}{A vector containing names of training datasets}
  \item{data.types}{A vector of molecular datatypes to load. Defaults to c('mRNA')}
  \item{subnets.file.flattened}{File containing all the binary ineractions derived from pathway-derived networks}
  \item{subset}{A list with a Field and Entry component specifying a subset of patients to be selected whose annotation Field matches Entry}
}
\value{Returns a list of matrices for each of the data types. Matrices contain nodes HR/P, edges HR and edges P.}
\author{Syed Haider & Paul C. Boutros}
\examples{
options("warn" = -1);
program.data <- get.program.defaults(networks.database = "test");
data.directory <- program.data[["test.data.dir"]];
subnets.file.flattened <- program.data[["subnets.file.flattened"]];
coef.nodes.edges <- calculate.network.coefficients(
  data.directory = data.directory,
  output.directory = ".",
  training.datasets = c("Breastdata1"),
  data.types = c("mRNA"),
  subnets.file.flattened = subnets.file.flattened
  );
}
\keyword{survival}

