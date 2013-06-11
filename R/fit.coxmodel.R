fit.coxmodel <- function(groups, survobj, stages = NA, rounding = 3, other.data = NULL) {

	# verify we have appropriate information
	if (nrow(survobj) != length(groups)) { return( rep(NA,5) ); }

	# fit a Cox model
	if (!is.null(other.data)) {
		# add all covariates
		other.data <- data.frame(groups, other.data);
		coxmodel <- coxph(survobj ~., data = other.data);
		}
	else {	
		if (nrow(survobj) == length(stages)) {
				coxmodel <- coxph(survobj ~ groups + stages);
				}
		else {
			coxmodel <- coxph(survobj ~ groups);
			}
		}

	# extract summary characteristics from the Cox fits (HR, 95% CIs, p)
	coxmodel  <- summary(coxmodel);
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

	# return results
	return( c(this.hr, this.95l, this.95u, this.pval, this.n) );

	}
