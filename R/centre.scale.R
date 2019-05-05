#' Centre and scale a data matrix
#' 
#' Centre and scale a data matrix. Scaling is done on each column separately
#' 
#' 
#' @param x A sample by feature data matrix
#' @param centre.data A character string specifying the centre value to be used for 
#' scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
#' e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
#' as well as for dichotomising data for estimating univariate betas from Cox model.
#' Defaults to 'median'
#' @return A centred and scaled data matrix
#' @author Syed Haider
#' @keywords center centre scale
#' @examples
#' 
#' tmp <- matrix(data = rnorm(100, 10, 2), nrow = 20);
#' tmp.scaled.median <- centre.scale.dataset(x = tmp);
#' tmp.scaled.mean <- centre.scale.dataset(x = tmp, centre.data = "mean");
#' tmp.scaled.custom <- centre.scale.dataset(x = tmp, centre.data = 0.3);
#' 
#' @export centre.scale.dataset
centre.scale.dataset <- function(x = NULL, centre.data = "median") {

	if (is.na(centre.data) || is.null(centre.data)) {
		stop("\ninvalid centre.data specified, please check documentation for valid values");
		}
	else if(centre.data == "median") {
		x.scaled <- apply(x, 2, FUN = function(y) { (y-median(y, na.rm = T))/stats::mad(y, na.rm = T) } );
		# for testing backwards compatibility
		# x.scaled <- scale(x);
		}
	else if (centre.data == "mean") {
		x.scaled <- scale(x);
		}
	else if (!is.na(as.numeric(centre.data))) {
		x.scaled <- scale(x, center = rep(as.numeric(centre.data), ncol(x)));
		}
	else {
		stop("\ninvalid centre.data specified, please check documentation for valid values");
		}

	return( x.scaled );
	}
