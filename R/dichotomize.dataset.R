#' Dichotomize a single dataset
#' 
#' Split a dataset into two groups by median-dichotomization
#' 
#' 
#' @param x A vector of values to be dichotomized
#' @param split.at An character string or a numeric value that is be used to dichotomize. 
#' Valid values are: 'median', 'mean', or a user defined numeric threshold. Defaults to
#' 'median'
#' @return A vector of the data dichotomized onto a 0/1 (low/high) scale.
#' @author Syed Haider & Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' tmp <- rnorm(100);
#' tmp.groups.median <- dichotomize.dataset(tmp);
#' tmp.groups.mean <- dichotomize.dataset(tmp, split.at = "mean");
#' tmp.groups.custom <- dichotomize.dataset(tmp, split.at = 0.3);
#' 
#' @export dichotomize.dataset
dichotomize.dataset <- function(x, split.at = "median") {

	if (is.na(split.at) || is.null(split.at)) {
		stop("\ninvalid split.at specified, please check documentation for valid values");
		}
	else if (split.at == "median") {
		split.at <- median(x, na.rm = TRUE);
		}
	else if (split.at == "mean") {
		split.at <- mean(x, na.rm = TRUE);
		}
	else if (!is.na(as.numeric(split.at))) {
		split.at <- as.numeric(split.at);
		}
	else {
		stop("\ninvalid split.at specified, please check documentation for valid values");
		}

	return( as.numeric( x > split.at ) );
	}
