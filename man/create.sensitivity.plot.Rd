% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.sensitivity.plot.R
\name{create.sensitivity.plot}
\alias{create.sensitivity.plot}
\title{Plots sensitivity analysis for class label dichotomization at supplied
survtime cutoffs}
\usage{
create.sensitivity.plot(
  riskscore = NULL,
  riskgroup = NULL,
  survtime = NULL,
  survstat = NULL,
  survtime.cutoffs = c(seq(5, 10, 1)),
  output.directory = ".",
  file.stem = NULL,
  main.title = "",
  resolution = 100
)
}
\arguments{
\item{riskscore}{A vector containing predicted risk scores}

\item{riskgroup}{A vector containing dichotomized risk groups}

\item{survtime}{A vector containing survival time of the samples}

\item{survstat}{A vector containing survival status of the samples}

\item{survtime.cutoffs}{A vector containing cutoff time points used to
dichotomize patients into low- and high-risk groups}

\item{output.directory}{Path to the output folder where intermediate and
results files will be saved}

\item{file.stem}{A string containing base name for image and text files
produced by this method}

\item{main.title}{A string specifying main title of the image}

\item{resolution}{A numeric value specifying resolution of the tiff image of
KM survival curves. Defaults to 100}
}
\value{
The sensitivity analysis plots are stored under
\code{output.directory}/graphs/. The sensitivity analysis results are stored
under \code{output.directory}/output/
}
\description{
A method to computer sensitivity, specificity and accuracy at all the
survtime cutoff steps provided
}
\author{
Syed Haider
}
\keyword{accuracy}
\keyword{sensitivity}
\keyword{specificity}
\keyword{survival}
