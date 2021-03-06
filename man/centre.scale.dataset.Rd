% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centre.scale.R
\name{centre.scale.dataset}
\alias{centre.scale.dataset}
\title{Centre and scale a data matrix}
\usage{
centre.scale.dataset(x = NULL, centre.data = "median")
}
\arguments{
\item{x}{A sample by feature data matrix}

\item{centre.data}{A character string specifying the centre value to be used for 
scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
as well as for dichotomising data for estimating univariate betas from Cox model.
Defaults to 'median'}
}
\value{
A centred and scaled data matrix
}
\description{
Centre and scale a data matrix. Scaling is done on each column separately
}
\examples{

tmp <- matrix(data = rnorm(100, 10, 2), nrow = 20);
tmp.scaled.median <- centre.scale.dataset(x = tmp);
tmp.scaled.mean <- centre.scale.dataset(x = tmp, centre.data = "mean");
tmp.scaled.custom <- centre.scale.dataset(x = tmp, centre.data = 0.3);

}
\author{
Syed Haider
}
\keyword{center}
\keyword{centre}
\keyword{scale}
