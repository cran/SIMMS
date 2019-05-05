#' Applies survdiff function
#' 
#' Applies survdiff on different prognoses groups and computes Logrank P using
#' chisquare statistics.
#' 
#' 
#' @param groups Grouping of patients (passed directly to survdiff, so factors
#' & continuous variables are okay)
#' @param survobj An object of class Surv (from the survival package) --
#' patient ordering needs to be identical as for groups
#' @return A vector containing: Chisq, degrees of freedom (DOF) and Logrank
#' P-value.
#' @author Syed Haider
#' @keywords survival
#' @examples
#' 
#' survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
#' survstat <- sample(c(0,1), 100, replace = TRUE);
#' survobj <- Surv(survtime, survstat);
#' groups <- sample(c('A','B'), 100, replace = TRUE);
#' get.chisq.stats(
#'   groups = as.factor(groups),
#'   survobj = survobj
#'   );
#' 
#' @export get.chisq.stats
get.chisq.stats <- function(groups, survobj) {

	# verify we have appropriate information
	if (nrow(survobj) != length(groups)) { return( rep(NA, 3) ); }

	# survdiff between the groups
	survdiff.obj <- survdiff(survobj ~ groups);

	# compute the logrank P
	logrank.p <- 1 - pchisq( survdiff.obj$chisq, df = (length(survdiff.obj$n) - 1) );

	# return results
	return(
		c(
			"Chisq" = survdiff.obj$chisq,
			"DOF" = (length(survdiff.obj$n) - 1),
			"LogRank.P" = logrank.p
			)
		);
	}
