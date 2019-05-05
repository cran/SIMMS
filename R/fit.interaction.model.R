#' Cox model two features separately and together
#' 
#' Using a meta-analysis dataset take two features and Cox model them
#' separately and together and extract HRs and p-values.
#' 
#' The interaction model compares cases where feature1 and feature2 concord
#' (both high or both low) to those where they do not. That is, the model is y
#' = x1 + x2 + (x1 == x2) and not the typical y = x1 + x2 + x1:x2
#' 
#' @param feature1 String indicate what feature (gene/probe/etc.) should be
#' extracted for analysis
#' @param feature2 String indicate what feature (gene/probe/etc.) should be
#' extracted for analysis
#' @param expression.data A list where each component is an expression matrix
#' (patients = columns, features = rows) for a different dataset
#' @param survival.data A list where each component is an object of class Surv
#' @param data.type.ordinal Logical indicating whether to treat this datatype
#' as ordinal. Defaults to FALSE
#' @param centre.data A character string specifying the centre value to be used for 
#' scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
#' e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
#' as well as for dichotomising data for estimating univariate betas from Cox model.
#' Defaults to 'median'
#' @return Returns a vector of six elements containing (HR,P) pairs for
#' feature1, feature2, and the interaction
#' @author Syed Haider & Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' data.dir <- get.program.defaults()[["test.data.dir"]];
#' data.types <- c("mRNA");
#' x1 <- load.cancer.datasets(
#'   datasets.to.load = c('Breastdata1'),
#'   data.types = data.types,
#'   data.directory = data.dir
#'   );
#' x2 <- fit.interaction.model(
#'   feature1 = "1000_at", 
#'   feature2 = "2549_at",
#'   expression.data = x1$all.data[[data.types[1]]],
#'   survival.data = x1$all.survobj
#'   );
#' 
#' @export fit.interaction.model
fit.interaction.model <- function(feature1, feature2, expression.data, survival.data, data.type.ordinal = FALSE, centre.data = "median") {

	groups1 <- SIMMS::dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature1,
		other.data = NULL,
		data.type.ordinal = data.type.ordinal,
		centre.data = centre.data
		);

	groups2 <- SIMMS::dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature2,
		other.data = NULL,
		data.type.ordinal = data.type.ordinal,
		centre.data = centre.data
		);

	# fit the interaction Cox model
	interaction.HR <- interaction.P <- NA;
	coxmodel <- NULL;
	survobj <- Surv(groups1$survtime, groups1$survstat);

	# handle the all-missing case smoothly
	if (all(is.na(groups1$groups * groups2$groups)) || 
		length(groups1$survtime) != length(groups1$groups) ||
		length(unique(groups1$groups[which(!is.na(groups1$groups))])) < 2 ||
		length(unique(groups2$groups[which(!is.na(groups2$groups))])) < 2
		) {
		cat('\nskipping interaction model, as all-missing from one or more feature groups or improperly formatted data (see below)');
		cat('\n\tunique levels in feature ', feature1, ":", length(unique(groups1$groups[which(!is.na(groups1$groups))])));
		cat('\n\tunique levels in feature ', feature2, ":", length(unique(groups2$groups[which(!is.na(groups2$groups))])));
		}
	else {

		# levels not specified as the beta of interaction terms remains unchanged regardless of 
		# which groups is used as base-line for groups1 and groups2
		tryCatch(
			expr = {
				coxmodel <- coxph(
					survobj ~ factor(groups1$groups) + factor(groups2$groups) + as.numeric(groups1$groups == groups2$groups)
					)
				},
			error = function(ex) {
				cat("\nInteraction model failed to converge (a known coxph issue) for features: ", feature1, " and ", feature2);
				}
			);

		# extract summary statistics from cox model
		coxmodel <- summary(coxmodel);

		# check fail to converge case
		if ("Class" %in% names(coxmodel) && "Mode" %in% names(coxmodel)
			&& coxmodel[["Class"]] == "NULL" && coxmodel[["Mode"]] == "NULL") {
			# no extra warning needed at this stage, but keep this placeholders for unforeseen coxph cases
			cat("");
			}
		else {
			interaction.HR <- as.numeric(coxmodel$coefficients["as.numeric(groups1$groups == groups2$groups)", 2]);
			interaction.P  <- as.numeric(coxmodel$coefficients["as.numeric(groups1$groups == groups2$groups)", 5]);
			}
		}

	# fit the two univariate models
	survival1 <- SIMMS::calculate.meta.survival(
		feature.name = feature1,
		expression.data = expression.data,
		survival.data = survival.data,
		data.type.ordinal = data.type.ordinal,
		centre.data = centre.data
		);

	survival2 <- SIMMS::calculate.meta.survival(
		feature.name = feature2,
		expression.data = expression.data,
		survival.data = survival.data,
		data.type.ordinal = data.type.ordinal,
		centre.data = centre.data
		);

	# return the survival statistics
	return(
		list(
			"cox.uv.1" = survival1,
			"cox.uv.2" = survival2,
			"cox.int" = c(interaction.HR, interaction.P)
			)
		);

	}
