#' Dichotomize a single dataset
#' 
#' Split a dataset into two groups by median-dichotomization
#' 
#' 
#' @param x A vector of values to be dichotomized
#' @param split.at An optional value that can be used to dichotomize instead of
#' median
#' @return A vector of the data dichotomized onto a 0/1 (low/high) scale.
#' @author Syed Haider & Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' tmp <- data.frame(y = rnorm(100));
#' tmp$x <- dichotomize.dataset(tmp$y);
#' 
#' @export dichotomize.dataset
dichotomize.dataset <- function(x, split.at = NA) {
	if (is.na(split.at)) { split.at = median(x, na.rm = TRUE); }
	return( as.numeric( x > split.at ) );
	}


