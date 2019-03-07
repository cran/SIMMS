#' Fit a Cox proportional hazards model
#' 
#' Fit a Cox model (possibly with some linear adjustments) and return key
#' statistics about the fit.
#' 
#' 
#' @param groups Grouping of patients (passed directly to coxph, so factors &
#' continuous variables are okay)
#' @param survobj An object of class Surv (from the survival package) --
#' patient ordering needs to be identical as for groups
#' @param stages DEPRECATED!  Use other.data instead.
#' @param rounding How many digits of precision should be returned?
#' @param other.data A data-frame (or matrix?) of variables to be controlled in
#' the Cox model. If null, no adjustment is done.  No interactions are fit.
#' @param data.type.ordinal Logical indicating whether to treat this datatype
#' as ordinal. Defaults to FALSE
#' @return A list containing two elements. \code{cox.stats} containing a vector
#' or matrix: HR, lower 95\% CI of HR, upper 95\% CI of HR, P-value (for
#' groups), number of samples (total with group assignments, although some may
#' not be included in fit for other reasons so this is an upper-limit).
#' \code{cox.obj} containing coxph model object
#' @author Syed Haider & Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
#' survstat <- sample(c(0,1), 100, replace = TRUE);
#' survobj <- Surv(survtime, survstat);
#' groups <- sample(c('A','B'), 100, replace = TRUE);
#' fit.coxmodel(
#'   groups = as.factor(groups),
#'   survobj = survobj
#'   );
#' 
#' @export fit.coxmodel
fit.coxmodel <- function(groups, survobj, stages = NA, rounding = 3, other.data = NULL, data.type.ordinal = FALSE) {

	# verify we have appropriate information
	to.abort <- FALSE;
	if (nrow(survobj) != length(groups)) { 
		to.abort <- TRUE;
		}

	# transform groups into factors if ordinal variable
	groups.unique <- unique(groups[which(!is.na(groups))]);

	if (length(groups.unique) < 2) {
		to.abort <- TRUE;
		cat("\nskipping cox model as less than 2 levels found in groups");
		}

	if (to.abort) {
		return( 
			list(
				"cox.stats" = rep(NA,5),
			 	"cox.obj" = NA
				)
			);
		}

	# process data for cox model
	if (data.type.ordinal) {

		# and set baseline to 0 (this is hardcoded at the moement)
		if (0 %in% groups) {
			groups.levels <- c(0, setdiff(groups.unique, 0));
			}
		else {
			groups.levels <- c(groups.unique);
			}

		groups <- factor(groups, levels = groups.levels);
		}

	# fit a Cox model
	if (!is.null(other.data)) {
		# add all covariates
		other.data <- data.frame(groups, other.data);
		coxfit <- coxph(survobj ~., data = other.data);
		}
	else {	
		if (nrow(survobj) == length(stages)) {
			coxfit <- coxph(survobj ~ groups + stages);
			}
		else {
			coxfit <- coxph(survobj ~ groups);
			}
		}

	# extract summary characteristics from the Cox fits (HR, 95% CIs, p)
	coxmodel  <- summary(coxfit);

	if (data.type.ordinal) {
		return.stats <- coxmodel$conf.int[, c(1, 3, 4), drop = FALSE]; # add HR, 95% CIs
		return.stats <- cbind(
			return.stats, 
			"p" = coxmodel$coef[, 5], # add p
			"n" = rep(length(groups[!is.na(groups)]), nrow(return.stats)) # add n
			);
		colnames(return.stats)[1:3] <- c("HR", "CI95L", "CI95U");
		rownames(return.stats) <- as.character(gsub("groups", "", rownames(return.stats)));

		# update n for each group
		return.stats[, "n"] <- table(groups)[rownames(return.stats)];
		}
	else {
		this.hr   <- coxmodel$conf.int[1,1];
		this.95l  <- coxmodel$conf.int[1,3];
		this.95u  <- coxmodel$conf.int[1,4];
		this.pval <- coxmodel$coef[1,5];
		this.n    <- length(groups[!is.na(groups)]);

		# round the major values to a few decimal places
		if (rounding) {
			this.hr  <- round(this.hr,  digits = rounding);
			this.95l <- round(this.95l, digits = rounding);
			this.95u <- round(this.95u, digits = rounding);
			}

		return.stats <- c(this.hr, this.95l, this.95u, this.pval, this.n);
		}

	# return results
	return(
		list(
			"cox.stats" = return.stats,
		 	"cox.obj" = coxfit
			)
		);

	}
