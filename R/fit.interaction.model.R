fit.interaction.model <- function(feature1, feature2, expression.data, survival.data) {

	groups1 <- dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature1,
		other.data = NULL
		);

	groups2 <- dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature2,
		other.data = NULL
		);

	# handle the all-missing case smoothly
	if (all(is.na(groups1$groups * groups2$groups)) || length(groups1$survtime) != length(groups1$groups)) {
		warning('All-missing from feature groups or improperly formatted data');
		return( rep(NA,6) );
		}

	# fit the interaction Cox model
	survobj <- Surv(groups1$survtime, groups1$survstat);
	coxmodel <- coxph(survobj ~ groups1$groups + groups2$groups + as.numeric(groups1$groups == groups2$groups));
	coxmodel <- summary(coxmodel);
	interaction.HR <- as.numeric(coxmodel$coefficients[3,2]);
	interaction.P  <- as.numeric(coxmodel$coefficients[3,5]);
	
	# fit the two univariate models
	survival1 <- calculate.meta.survival(
		feature.name = feature1,
		expression.data = expression.data,
		survival.data = survival.data
		);

	survival2 <- calculate.meta.survival(
		feature.name = feature2,
		expression.data = expression.data,
		survival.data = survival.data
		);
			
	# return the survival statistics
	return(
		c(
			survival1[c(1,4)],
			survival2[c(1,4)],
			interaction.HR,
			interaction.P
			)
		);
	
	}
