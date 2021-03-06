% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.survobj.R
\name{create.survobj}
\alias{create.survobj}
\title{Utility function for loading meta-analysis lists}
\usage{
create.survobj(annotation = NULL, truncate.survival = 100)
}
\arguments{
\item{annotation}{A patient annotation matrix (patients = rows) with (at
least) columns for survtime, survstat, and survtime.unit}

\item{truncate.survival}{A numeric value specifying survival truncation in
years. Defaults to 100 years which effectively means no truncation}
}
\value{
Returns an object of class Surv
}
\description{
Create Surv objects from an annotation-matrix with handling for different
time units.
}
\examples{

annotation.file <- paste(
  get.program.defaults()[["test.data.dir"]],
  "/Breastdata2/patient_annotation.txt", sep = ""
  );
annotation <- read.table(
  annotation.file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
  );

# select the appropriate survtime and survstat variable for this dataset
annotation$survstat      <- annotation[,'e.dfs'];
annotation$survtime      <- annotation[,'t.dfs'];
annotation$survtime.unit <- annotation[,'t.dfs.unit'];

# only keep samples with survival data
annotation <- annotation[!is.na(annotation$survstat) & !is.na(annotation$survstat),];

surv.obj <- create.survobj(annotation = annotation);

}
\author{
Paul C. Boutros
}
\keyword{survival}
