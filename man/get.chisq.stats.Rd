% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.chisq.stats.R
\name{get.chisq.stats}
\alias{get.chisq.stats}
\title{Applies survdiff function}
\usage{
get.chisq.stats(groups, survobj)
}
\arguments{
\item{groups}{Grouping of patients (passed directly to survdiff, so factors
& continuous variables are okay)}

\item{survobj}{An object of class Surv (from the survival package) --
patient ordering needs to be identical as for groups}
}
\value{
A vector containing: Chisq, degrees of freedom (DOF) and Logrank
P-value.
}
\description{
Applies survdiff on different prognoses groups and computes Logrank P using
chisquare statistics.
}
\examples{

survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
survstat <- sample(c(0,1), 100, replace = TRUE);
survobj <- Surv(survtime, survstat);
groups <- sample(c('A','B'), 100, replace = TRUE);
get.chisq.stats(
  groups = as.factor(groups),
  survobj = survobj
  );

}
\author{
Syed Haider
}
\keyword{survival}
