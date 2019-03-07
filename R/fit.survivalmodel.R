#' Trains a multivariate survival model
#' 
#' Trains a multivariate survival model and conducts feature selection using
#' both backward elimination and forward selection, independently. TO BE
#' DEPRECATED AND HAS BEEN REPLACED BY \code{create.classifier.multivariate}
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{feature.selection.datasets}, \code{training.datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files will be saved
#' @param feature.selection.datasets A vector containing names of datasets used
#' for feature selection in function \code{derive.network.features()}
#' @param feature.selection.p.threshold One of the P values that were used for
#' feature selection in function \code{derive.network.features()}. This
#' function does not support vector of P values as used in
#' \code{derive.network.features()} for performance reasons
#' @param training.datasets A vector containing names of training datasets to
#' be used to train multivariate survival model
#' @param top.n.features A numeric value specifying how many top ranked
#' features will be used to train the multivariate survival model
#' @param models A character vector specifying which models ('1' = N+E, '2' =
#' N, '3' = E) to run
#' @return The output files are stored under \code{output.directory}/output/
#' @author Syed Haider
#' @seealso \code{create.classifier.multivariate}
#' @keywords survival
#' @examples
#' 
#' # see package's main documentation
#' 
#' @export fit.survivalmodel
fit.survivalmodel <- function(data.directory = ".", output.directory = ".", feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05, training.datasets = NULL, top.n.features = 25, models = c("1", "2", "3")) {

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# program data files and initialise variables
	all.feature.selection.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.training.names <- paste(sort(training.datasets), collapse="_");
	top.subnets <- list();
	all.training.data <- NULL;
	stepAIC.results <- NULL;

	# lets read top subnets for the give feature.selection.datasets
	for (model in models) {
		top.subnets[[model]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/top_subnets_score__TRAINING_", all.feature.selection.names, "__model_", model, "__PV_", feature.selection.p.threshold, ".txt", sep = ""),
				header = TRUE,
				row.names = 1,
				sep = "\t"
				)
			);

		# join the training datasets together to make a large training cohort
		for (dataset in training.datasets) {
			all.training.data[[model]] <- cbind(
				all.training.data[[model]],
				as.matrix(
					read.table(
						file = paste(out.dir, "/patientwise_subnets_score__", dataset, "__TRAINING_", all.feature.selection.names, "__model_", model, ".txt", sep = ""),
						header = TRUE,
						row.names = 1,
						sep = "\t"
						)
					)
				);
			}

		}

	# run stepAIC for feature subset selection (backward and forward)
	for (model in models) {

		# in case requested number of features are greater than available features
		if (top.n.features > nrow(top.subnets[[model]])) {
			stop("\ntop.n.features larger than available features i-e ", nrow(top.subnets[[model]]));
			}

		explanatory.variables <- rownames(top.subnets[[model]])[1:top.n.features];
		model.formula <- list();
		all.training.data.trans <- t(all.training.data[[model]]);
		model.formula[["backward"]] <- as.formula(
			paste(
				"Surv(survtime, survstat) ~ ",
				paste(
					explanatory.variables,
					collapse = "+"
					)
				)
			);
		model.formula[["forward"]] <- as.formula(
			paste(
				"Surv(survtime, survstat) ~ ",
				paste(
					" 1",
					collapse = ""
					)
				)
			);

		for (direction in c("forward", "backward")) {
			cat("\n*****************  starting fit\tModel ", model, direction, "  *****************\n");
			print(model.formula[[direction]]);
			stepAIC.results <- stepAIC(
				coxph(model.formula[[direction]], data = as.list(as.data.frame(all.training.data.trans))),
				direction = direction,
				trace = TRUE,
				scope = list(
					lower = coxph(model.formula[["forward"]], data = as.list(as.data.frame(all.training.data.trans))),
					upper = coxph(model.formula[["backward"]], data = as.list(as.data.frame(all.training.data.trans))) 
					)
				);

			# lets store the results of each model (summary) to file
			model.summary <- summary(stepAIC.results);
			if (is.null(model.summary$conf.int)) {
				stop("\nNull Model, no variables selected after stepAIC");
				}
			res.matrix <- matrix(
				data = NA,
				ncol = 5,
				nrow = nrow(model.summary$conf.int),
				dimnames = list(
					rownames(model.summary$conf.int),
					c("HR", "95l", "95u", "p-val", "beta")
					)
				);

			# lets extract the feature name (i.e. subnetwork name)
			for(feature.name in rownames(model.summary$conf.int)) {
				#store HR, 95l, 95u, p-val
				res.matrix[feature.name, 1] <- model.summary$conf.int[feature.name, 1];
				res.matrix[feature.name, 2] <- model.summary$conf.int[feature.name, 3];
				res.matrix[feature.name, 3] <- model.summary$conf.int[feature.name, 4];
				res.matrix[feature.name, 4] <- model.summary$coef[feature.name, 5];
				res.matrix[feature.name, 5] <- model.summary$coef[feature.name, 1];
				}

			# lets store it to the file system
			cat("\nstoring StepAIC results for: ", model);
			write.table(
				res.matrix,
				file = paste(out.dir, "/beta__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			cat("\nstoring AIC TRACE, the last one being the AIC of the final model");
			write.table(
				as.matrix(stepAIC.results$anova[1:6]),
				file = paste(out.dir, "/AIC__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);
			}
		}

	# return(all.training.data);
	# PCB: why comment this out? [its only for debugging :)]
	}
