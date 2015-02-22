\name{prepare.training.validation.datasets}
\alias{prepare.training.validation.datasets}
\title{Prepare training and validation datasets}
\description{Computes per-patient pathway-derived network impact scores across all input datasets, independently}
\usage{
prepare.training.validation.datasets(
  data.directory = ".", output.directory = ".",
  data.types = c("mRNA"), data.types.ordinal = c("cnv"),
  min.ordinal.threshold = c("cnv" = 3), p.threshold = 1, 
  feature.selection.datasets = NULL, datasets = NULL, 
  truncate.survival = 100, networks.database = "default", 
  write.normed.datasets = TRUE, subset = NULL
  );
}
\arguments{
  \item{data.directory}{Path to the directory containing datasets as specified by \code{datasets}}
  \item{output.directory}{Path to the output folder where intermediate and results files will be saved}
  \item{data.types}{A vector of molecular datatypes to load. Defaults to c('mRNA')}
  \item{data.types.ordinal}{A vector of molecular datatypes to be treated as ordinal. Defaults to c('cnv')}
  \item{min.ordinal.threshold}{A named vector specifying minimum percent threshold for each ordinal data type to be used prior to estimating coefficients. Coefficient for features not satisfying minimum threshold will not be estimated, and set to 0. Defaults to cnv threshold as 3 percent}
  \item{p.threshold}{P value threshold to be applied for selecting univariate prognostic features. Defaults to 1}
  \item{feature.selection.datasets}{A vector containing names of datasets used for feature selection in function \code{derive.network.features()}}
  \item{datasets}{A vector containing names of all the datasets to be later used for training and validation purposes}
  \item{truncate.survival}{A numeric value specifying survival truncation in years. Defaults to 100 years which effectively means no truncation}
  \item{networks.database}{Name of the pathway networks database. Default to NCI PID/Reactome/Biocarta i-e "default"}
  \item{write.normed.datasets}{A toggle to control whether processed mRNA and survival data should be written to file}
  \item{subset}{A list with a Field and Entry component specifying a subset of patients to be selected whose annotation Field matches Entry}
}
\value{The output files are stored under \code{output.directory}/output/}
\author{Syed Haider}
\examples{
# get data directory 
data.directory <- get.program.defaults()[["test.data.dir"]];

# initialise params
output.directory <- ".";
data.types <- c("mRNA");
feature.selection.datasets <- c("Breastdata1");
training.datasets <- c("Breastdata1");
validation.datasets <- c("Breastdata1", "Breastdata2");

# preparing training and validation datasets.
# Normalisation & patientwise subnet feature scores
prepare.training.validation.datasets(
  data.directory = data.directory,
  output.directory = output.directory,
  data.types =  data.types,
  feature.selection.datasets = feature.selection.datasets,
  datasets = unique(c(training.datasets, validation.datasets)),
  networks.database = "test"
  );
}
\keyword{IO}
