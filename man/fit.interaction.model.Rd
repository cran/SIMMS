\name{fit.interaction.model}
\alias{fit.interaction.model}
\title{Cox model two features separately and together}
\description{Using a meta-analysis dataset take two features and Cox model them separately and together and extract HRs and p-values.}
\usage{
fit.interaction.model(
  feature1, feature2,
  expression.data, survival.data,
  data.type.ordinal = FALSE
  );
}
\arguments{
  \item{feature1}{String indicate what feature (gene/probe/etc.) should be extracted for analysis}
  \item{feature2}{String indicate what feature (gene/probe/etc.) should be extracted for analysis}
  \item{expression.data}{A list where each component is an expression matrix (patients = columns, features = rows) for a different dataset}
  \item{survival.data}{A list where each component is an object of class Surv}
  \item{data.type.ordinal}{Logical indicating whether to treat this datatype as ordinal. Defaults to FALSE}
}
\details{
The interaction model compares cases where feature1 and feature2 concord (both high or both low) to those where they do not.
That is, the model is y = x1 + x2 + (x1 == x2) and not the typical y = x1 + x2 + x1:x2
}
\value{Returns a vector of six elements containing (HR,P) pairs for feature1, feature2, and the interaction}
\author{Paul C. Boutros}
\examples{
data.dir <- get.program.defaults()[["test.data.dir"]];
data.types <- c("mRNA");
x1 <- load.cancer.datasets(
  datasets.to.load = c('Breastdata1'),
  data.types = data.types,
  data.directory = data.dir
  );
x2 <- fit.interaction.model(
  feature1 = "1000_at", 
  feature2 = "2549_at",
  expression.data = x1$all.data[[data.types[1]]],
  survival.data = x1$all.survobj
  );
}
\keyword{survival}
