\name{create.classifier.multivariate}
\alias{create.classifier.multivariate}
\title{Trains and tests a multivariate survival model}
\description{Trains a model on training datasets. Predicts the risk score for all the training & datasets, independently. This function also predicts the risk score for combined training datasets cohort and validation datasets cohort. The risk score estimation is done by multivariate models fit by \code{fit.survivalmodel}. The function also predicts risk scores for each of the \code{top.n.features} independently.}
\usage{
create.classifier.multivariate(
  data.directory = ".", output.directory = ".",
  feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05,
  training.datasets = NULL, validation.datasets = NULL,
  top.n.features = 25, models = c("1", "2", "3"),
  learning.algorithms = c("backward", "forward"),
  alpha.glm = c(1), k.fold.glm = 10, seed.cv.glm = 51214,
  cores.glm = 1
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
  \item{learning.algorithms}{A character vector specifying which learning algorithm to be used for model fitting and feature selection. Defaults to c('backward', 'forward'). Available options are: c('backward', 'forward', 'glm')}
  \item{alpha.glm}{A numeric vector specifying elastic-net mixing parameter alpha, with range alpha raning from [0,1]. 1 for LASSO (default) and 0 for ridge. For multiple values of alpha, most optimal value is selected through cross validation on training set}
  \item{k.fold.glm}{A numeric value specifying k-fold cross validation if glm was chosen in \code{learning.algorithms}}
  \item{seed.cv.glm}{A numeric value specifying seed for k-fold cross validation if glm was chosen in \code{learning.algorithms}}
  \item{cores.glm}{An integer value specifying number of cores to be used for glm if it was chosen in \code{learning.algorithms}}
}
\value{The output files are stored under \code{output.directory}/output/}
\author{Syed Haider & Vincent Stimper}
\examples{
# see package's main documentation
}
\keyword{survival}
