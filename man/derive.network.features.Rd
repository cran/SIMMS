\name{derive.network.features}
\alias{derive.network.features}
\title{Derive univariate features from pathway-derived networks}
\description{This function fits Cox model to features as well as interaction between features. The coefficients of features are subsequently used to compute impact score of each of the pathway-derived networks.}
\usage{
derive.network.features(
  data.directory = ".", output.directory = ".", 
  data.types = c("mRNA"), data.types.ordinal = c("cnv"),
  feature.selection.fun = "calculate.network.coefficients",
  feature.selection.datasets = NULL, 
  feature.selection.p.thresholds = c(0.05), truncate.survival = 100,
  networks.database = "default", subset = NULL, ...
  );
}
\arguments{
  \item{data.directory}{Path to the directory containing datasets as specified by \code{feature.selection.datasets}}
  \item{output.directory}{Path to the output folder where intermediate and results files will be saved}
  \item{data.types}{A vector of molecular datatypes to load. Defaults to c('mRNA')}
  \item{data.types.ordinal}{A vector of molecular datatypes to be treated as ordinal. Defaults to c('cnv')}
  \item{feature.selection.fun}{Name of the function to be used to estimate network coefficients. Defaults to 'calculate.network.coefficients'}
  \item{feature.selection.datasets}{A vector containing names of training datasets to be used to compute cox statistics}
  \item{feature.selection.p.thresholds}{A vector containing P values to be used as threshold for including features into overall impact score of a network}
  \item{truncate.survival}{A numeric value specifying survival truncation in years. Defaults to 100 years which effectively means no truncation}
  \item{networks.database}{Name of the pathway networks database. Default to NCI PID/Reactome/Biocarta i-e "default"}
  \item{subset}{A list with a Field and Entry component specifying a subset of patients to be selected from each dataset whose annotation Field matches Entry}
  \item{...}{other params to be passed on to user-defined method for estimating coefficients of network features}
}
\value{The output files are stored under \code{data.directory}/output/}
\author{Syed Haider}
\examples{
options("warn" = -1);

# get data directory 
data.directory <- get.program.defaults(networks.database = "test")[["test.data.dir"]];

# initialise params
output.directory <- ".";
data.types <- c("mRNA");
feature.selection.datasets <- c("Breastdata1");
feature.selection.p.thresholds <- c(0.05);

# estimate network coefficients for all the subnet features
derive.network.features(
  data.directory = data.directory,
  output.directory = output.directory,
  data.types = data.types,
  feature.selection.fun = "calculate.network.coefficients",
  feature.selection.datasets = feature.selection.datasets,
  feature.selection.p.thresholds = feature.selection.p.thresholds,
  networks.database = "test"
  );
}
\keyword{FeatureSelection}
