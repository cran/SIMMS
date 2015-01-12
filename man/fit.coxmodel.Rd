\name{fit.coxmodel}
\alias{fit.coxmodel}
\title{Fit a Cox proportional hazards model}
\description{Fit a Cox model (possibly with some linear adjustments) and return key statistics about the fit.}
\usage{
fit.coxmodel(
  groups, survobj, stages = NA,
  rounding = 3, other.data = NULL,
  data.type.ordinal = FALSE
  );
}
\arguments{
  \item{groups}{Grouping of patients (passed directly to coxph, so factors & continuous variables are okay)}
  \item{survobj}{An object of class Surv (from the survival package) -- patient ordering needs to be identical as for groups}
  \item{stages}{DEPRECATED!  Use other.data instead.}
  \item{rounding}{How many digits of precision should be returned?}
  \item{other.data}{A data-frame (or matrix?) of variables to be controlled in the Cox model. If null, no adjustment is done.  No interactions are fit.}
  \item{data.type.ordinal}{Logical indicating whether to treat this datatype as ordinal. Defaults to FALSE}
}
\value{A list containing two elements. \code{cox.stats} containing a vector or matrix: HR, lower 95\% CI of HR, upper 95\% CI of HR, P-value (for groups), number of samples (total with group assignments, although some may not be included in fit for other reasons so this is an upper-limit). \code{cox.obj} containing coxph model object}
\author{Paul C. Boutros}
\examples{
survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
survstat <- sample(c(0,1), 100, replace = TRUE);
survobj <- Surv(survtime, survstat);
groups <- sample(c('A','B'), 100, replace = TRUE);
fit.coxmodel(
  groups = as.factor(groups),
  survobj = survobj
  );
}
\keyword{survival}
