\name{create.KM.plot}
\alias{create.KM.plot}
\title{Plots Kaplan-meier survival curve for a given risk grouping & survival params}
\description{A generic method to plot KM curves}
\usage{
create.KM.plot(
  riskgroup = NULL, survtime =  NULL, survstat = NULL,
  file.name = NULL, main.title = "", resolution = 100
  );
}
\arguments{
  \item{riskgroup}{A vector containing dichotomized risk groups}
  \item{survtime}{A vector containing survival time of the samples}
  \item{survstat}{A vector containing survival status of the samples}
  \item{file.name}{A string containing full qualified path of the output tiff file}
  \item{main.title}{A string specifying main title of the image}
  \item{resolution}{A numeric value specifying resolution of the tiff image of KM survival curves. Defaults to 100}
}
\value{The KM survival curves are stored under \code{output.dir}/graphs/}
\author{Syed Haider}
\keyword{survival,Kaplan-meier}
