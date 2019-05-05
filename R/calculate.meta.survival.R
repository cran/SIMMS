#' Fit a meta-analytic Cox proportional hazards model to a single feature
#' 
#' Takes a meta-analysis data object and fits a Cox proportional hazards model
#' (possibly with adjustment for some specific covariates) by
#' median-dichotomizing patients within each individual dataset.
#' 
#' 
#' @param feature.name Character indicate what feature (gene/probe/etc.) should
#' be extracted for analysis
#' @param expression.data A list where each component is an expression matrix
#' (patients = columns, genes = rows) for a different dataset
#' @param survival.data A list where each component is an object of class Surv
#' @param rounding How many digits after the decimal place to include
#' @param other.data A list of other covariates to be passed to the Cox model
#' (all elements in this list are used
#' @param data.type.ordinal Logical indicating whether to treat this datatype
#' as ordinal. Defaults to FALSE
#' @param centre.data A character string specifying the centre value to be used for 
#' scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
#' e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
#' as well as for dichotomising data for estimating univariate betas from Cox model.
#' Defaults to 'median'
#' @return Returns a vector containing the HR, p-value, n, and 95\% confidence
#' limits of the HR (see fit.coxmodel() for details)
#' @author Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' data.directory <- get.program.defaults()[["test.data.dir"]];
#' data.types <- c("mRNA");
#' x1 <- load.cancer.datasets(
#'   datasets.to.load = c('Breastdata1'),
#'   data.types = data.types,
#'   data.directory = data.directory
#'   );
#' x2 <- calculate.meta.survival(
#'   feature.name = "1000_at",
#'   expression.data = x1$all.data[[data.types[1]]],
#'   survival.data = x1$all.survobj
#'   );
#' 
#' @export calculate.meta.survival
calculate.meta.survival <- function(feature.name, expression.data, survival.data, rounding = 3, other.data = NULL, data.type.ordinal = FALSE, centre.data = "median") {

	# verify that we got appropriate input data
	to.abort <- FALSE;
	if ("list" != class(expression.data)) { to.abort <- TRUE; }
	if ("list" != class(survival.data)) { to.abort <- TRUE; }
	if (length(expression.data) != length(survival.data)) { to.abort <- TRUE; }

	# stop processing if we have bad data
	if (to.abort) {
		warning('data failed sanity-checking');
		return(
			list(
				"cox.stats" = rep(NA,5),
			 	"cox.obj" = NA
				)
			);
		}

	# dichotomize meta
	dichotomized.data <- SIMMS::dichotomize.meta.dataset(
		feature.name = feature.name,
		expression.data = expression.data,
		survival.data = survival.data,
		other.data = other.data,
		data.type.ordinal = data.type.ordinal,
		centre.data = centre.data
		);

	# handle all-NA values (i.e. feature not in the dataset)
	if (0 == length(dichotomized.data$groups) || 
		0 == length(dichotomized.data$survtime) || 
		0 == length(dichotomized.data$survstat) || 
		all(is.na(dichotomized.data$groups))
		) {
		cat('\n\tfeature not in the dataset: ', feature.name);
		return( 
			list(
				"cox.stats" = rep(NA,5),
			 	"cox.obj" = NA
				)
			);
		}

	# fit coxph model
	return(
		SIMMS::fit.coxmodel(
			groups = dichotomized.data$groups,
			survobj = Surv(dichotomized.data$survtime, dichotomized.data$survstat),
			rounding = rounding,
			other.data = other.data,
			data.type.ordinal = data.type.ordinal
			)
		);
	}
