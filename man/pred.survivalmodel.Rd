\name{pred.survivalmodel}
\alias{pred.survivalmodel}
\title{Trains a multivariate survival model}
\description{Predicts the risk score for all the training & datasets, independently. This function predicts the risk score for combined training datasets cohort and validation datasets cohort. The risk score computation is done for multivariate models fit by \code{fit.survivalmodel}. The function also predicts risk scores for each of the \code{top.n.features} independently.}
\usage{
pred.survivalmodel(
  data.directory = ".", output.directory = ".",
  feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05,
  training.datasets = NULL, validation.datasets = NULL,
  top.n.features = 25,  models = c("1", "2", "3"),
  truncate.survival = 100, write.risk.data = TRUE
  );
}
\arguments{
  \item{data.directory}{Path to the directory containing datasets as specified by \code{feature.selection.datasets}, \code{training.datasets}, \code{validation.datasets}}
  \item{output.directory}{Path to the output folder where intermediate and results files will be saved}
  \item{feature.selection.datasets}{A vector containing names of datasets used for feature selection in function \code{derive.network.features()}}
  \item{feature.selection.p.threshold}{One of the P values that were used for feature selection in function \code{derive.network.features()}. This function does not support vector of P values as used in \code{derive.network.features()} for performance reasons}
  \item{training.datasets}{A vector containing names of training datasets}
  \item{validation.datasets}{A vector containing names of validation datasets}
  \item{top.n.features}{A numeric value specifying how many top ranked features will be used for univariate survival modelling}
  \item{models}{A character vector specifying which of the models ('1' = N+E, '2' = N, '3' = E) to run}
  \item{truncate.survival}{A numeric value specifying survival truncation in years. Defaults to 100 years which effectively means no truncation}
  \item{write.risk.data}{A toggle to control whether risk scores and patient risk groups should be written to file}
}
\value{The output files are stored under \code{output.directory}/output/}
\author{Syed Haider}
\examples{
# see package's main documentation
}
\keyword{survival}
