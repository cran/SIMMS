\name{fit.survivalmodel}
\alias{fit.survivalmodel}
\title{Trains a multivariate survival model}
\description{Trains a multivariate survival model and conducts feature selection using both backward elimination and forward selection, independently. TO BE DEPRECATED AND HAS BEEN REPLACED BY \code{create.classifier.multivariate}}
\usage{
fit.survivalmodel(
  data.directory = ".", output.directory = ".",
  feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05,
  training.datasets = NULL, top.n.features = 25, models = c("1", "2", "3")
  );
}
\arguments{
  \item{data.directory}{Path to the directory containing datasets as specified by \code{feature.selection.datasets}, \code{training.datasets}}
  \item{output.directory}{Path to the output folder where intermediate and results files will be saved}
  \item{feature.selection.datasets}{A vector containing names of datasets used for feature selection in function \code{derive.network.features()}}
  \item{feature.selection.p.threshold}{One of the P values that were used for feature selection in function \code{derive.network.features()}. This function does not support vector of P values as used in \code{derive.network.features()} for performance reasons}
  \item{training.datasets}{A vector containing names of training datasets to be used to train multivariate survival model}
  \item{top.n.features}{A numeric value specifying how many top ranked features will be used to train the multivariate survival model}
  \item{models}{A character vector specifying which models ('1' = N+E, '2' = N, '3' = E) to run}
}
\value{The output files are stored under \code{output.directory}/output/}
\author{Syed Haider}
\seealso{
\code{create.classifier.multivariate}
}
\examples{
# see package's main documentation
}
\keyword{survival}
