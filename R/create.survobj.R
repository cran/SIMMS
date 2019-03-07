#' Utility function for loading meta-analysis lists
#' 
#' Create Surv objects from an annotation-matrix with handling for different
#' time units.
#' 
#' 
#' @param annotation A patient annotation matrix (patients = rows) with (at
#' least) columns for survtime, survstat, and survtime.unit
#' @param truncate.survival A numeric value specifying survival truncation in
#' years. Defaults to 100 years which effectively means no truncation
#' @return Returns an object of class Surv
#' @author Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' annotation.file <- paste(
#'   get.program.defaults()[["test.data.dir"]],
#'   "/Breastdata2/patient_annotation.txt", sep = ""
#'   );
#' annotation <- read.table(
#'   annotation.file,
#'   header = TRUE,
#'   row.names = 1,
#'   sep = "\t"
#'   );
#' 
#' # select the appropriate survtime and survstat variable for this dataset
#' annotation$survstat      <- annotation[,'e.dfs'];
#' annotation$survtime      <- annotation[,'t.dfs'];
#' annotation$survtime.unit <- annotation[,'t.dfs.unit'];
#' 
#' # only keep samples with survival data
#' annotation <- annotation[!is.na(annotation$survstat) & !is.na(annotation$survstat),];
#' 
#' surv.obj <- create.survobj(annotation = annotation);
#' 
#' @export create.survobj
create.survobj <- function(annotation = NULL, truncate.survival = 100) {

	# localize the survival data
	survtime <- annotation$survtime;
	survstat <- annotation$survstat;

	# check the annotation for each patient
	if (length(annotation$survtime.unit) == 0 | length(annotation$survtime) == 0) {
		warning('patient annotation missing');
		return(NA);
		}
	
	for (j in 1:nrow(annotation)) {
		if (is.na(annotation$survtime.unit[j]) | is.na(survtime[j])) { survtime[j] <- NA; }
		else if (annotation$survtime.unit[j] == "days")   { survtime[j] <- survtime[j] / 365.25; }
		else if (annotation$survtime.unit[j] == "weeks")  { survtime[j] <- survtime[j] / 52.18; }
		else if (annotation$survtime.unit[j] == "months") { survtime[j] <- survtime[j] / 12; }
		else if (annotation$survtime.unit[j] == "years")  { } # do nothing, this is the default}
		}

	# create the Surv object
	if (!is.numeric(survtime) || !is.numeric(survstat)) {
		survobj <- rep(NA, length(survtime));
		}
	else {
		# truncate survival
		survstat[survtime > truncate.survival] <- 0;
		survtime[survtime > truncate.survival] <- truncate.survival;

		# create survival object
		survobj <- Surv(survtime, survstat);
		}

	return(survobj);

	}
