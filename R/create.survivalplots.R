#' Plots Kaplan-meier survival curves
#' 
#' Plots Kaplan-meier survival curves for all the training & datasets,
#' independently as well as combined training datasets cohort and validation
#' datasets cohort. The function also plots KM survival curves for each of the
#' \code{top.n.features} independently.
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{training.datasets}, \code{validation.datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files were saved
#' @param training.datasets A vector containing names of training datasets
#' @param validation.datasets A vector containing names of validation datasets
#' @param top.n.features A numeric value specifying how many top ranked
#' features will be used for univariate survival modelling
#' @param learning.algorithms A character vector specifying which learning
#' algorithm to be used for model fitting and feature selection. Defaults to
#' c('backward', 'forward'). Available options are: c('backward', 'forward',
#' 'glm')
#' @param truncate.survival A numeric value specifying survival truncation in
#' years. Defaults to 100 years which effectively means no truncation
#' @param survtime.cutoffs A vector containing survival cutoff time points to
#' be used for dichotomization of patients into risk groups for senstivity
#' analysis
#' @param main.title A logical to specify plot's main title. Defaults to FASLE
#' @param KM.plotting.fun A string containing the name of the method to use for
#' plotting KM curves. Defaults to \code{create.KM.plot}
#' @param plot.univariate.data Logical to indicate whether to plot univariate
#' results for all subnetworks. Default to FALSE
#' @param plot.multivariate.data Logical to indicate whether to plot
#' multivariate results for all subnetworks. Defaults to TRUE
#' @param resolution A numeric value specifying resolution of the png images of
#' KM survival curves. Defaults to 100
#' @return The KM survival curves are stored under
#' \code{output.directory}/graphs/
#' @author Syed Haider
#' @keywords survival,Kaplan-meier
#' @examples
#' 
#' # see package's main documentation
#' 
#' @export create.survivalplots
create.survivalplots <- function(data.directory = ".", output.directory = ".", training.datasets = NULL, validation.datasets = NULL, top.n.features = 25, learning.algorithms = c("backward", "forward"), truncate.survival = 100, survtime.cutoffs = c(seq(5,10,1)), main.title = FALSE, KM.plotting.fun = "create.KM.plot", plot.univariate.data = FALSE, plot.multivariate.data = TRUE, resolution = 100) {

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# program data files and initialise variables
	all.training.names <- paste(sort(training.datasets), collapse="_");
	all.validation.names <- paste(sort(validation.datasets), collapse="_");
	models <- c("1", "2", "3");
	model.names <- c("N+E", "N", "E");
	all.riskgroups.data <- list();
	all.riskscores.data <- list();

	# lets read in the training & validation dataset riskgroups (multivariate models)
	if (plot.multivariate.data) {
		for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
			for (model in models) {
				for (learning.algorithm in learning.algorithms) {

					riskgroups.file <- paste(out.dir, "/riskgroups__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep="");

					if (file.exists(riskgroups.file)) {
						all.riskgroups.data[[dataset]][[model]][[learning.algorithm]] <- as.matrix(
							read.table(
								file = riskgroups.file,
								header = TRUE,
								row.names = 1,
								sep = "\t"
								)
							);
						}

					riskscores.file <- paste(out.dir, "/riskscores__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep="");

					if (file.exists(riskscores.file)) {
						all.riskscores.data[[dataset]][[model]][[learning.algorithm]] <- as.matrix(
							read.table(
								file = riskscores.file,
								header = TRUE,
								row.names = 1,
								sep = "\t"
								)
							);
						}
					}
				}
			}

		# lets do KM plotting (multivariate models)
		for (dataset in names(all.riskgroups.data)) {
			for (model in names(all.riskgroups.data[[dataset]])) {
				for (learning.algorithm in names(all.riskgroups.data[[dataset]][[model]])) {

					# call appropriate KM plotting function, user defined OR this package's default
					eval(
						call(
							KM.plotting.fun,
							riskgroup = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "riskgroup"],
							survtime = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "survtime"],
							survstat = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "survstat"],
							file.name = paste(graphs.dir, "KM__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, "__truncated_", truncate.survival, ".png", sep=""),
							main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ", ", learning.algorithm, ")", sep = ""), ""),
							resolution = resolution
							)
						);

					# call sensitivity analysis function
					create.sensitivity.plot(
						riskscore = all.riskscores.data[[dataset]][[model]][[learning.algorithm]][, "riskscore"],
						riskgroup = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "riskgroup"],
						survtime = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "survtime"],
						survstat = all.riskgroups.data[[dataset]][[model]][[learning.algorithm]][, "survstat"],
						survtime.cutoffs = survtime.cutoffs,
						output.directory = output.directory,
						file.stem = paste("sensitivity_analysis__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, sep=""),
						main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ", ", learning.algorithm, ")", sep = ""), ""),
						resolution = resolution
						);
					}
				}
			}
		}

	all.riskgroups.data <- list();
	all.riskscores.data <- list();

	# lets read in the training & validation dataset riskgroups (univariate grouping)
	if (plot.univariate.data) {
		for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
			for (model in models) {

				riskgroups.file <- paste(out.dir, "riskgroups_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep="");
				if (file.exists(riskgroups.file)) {
					all.riskgroups.data[[dataset]][[model]] <- as.matrix(
						read.table(
							file = riskgroups.file,
							header = TRUE,
							row.names = 1,
							sep = "\t"
							)
						);
					}

				riskscores.file <- paste(out.dir, "riskscores_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep="");
				if (file.exists(riskscores.file)) {
					all.riskscores.data[[dataset]][[model]] <- as.matrix(
						read.table(
							file = riskscores.file,
							header = TRUE,
							row.names = 1,
							sep = "\t"
							)
						);
					}
				}
			}

		# let's do KM plotting (univariate grouping)
		for (dataset in names(all.riskgroups.data)) {
			for (model in names(all.riskgroups.data[[dataset]])) {
				feature.names <- colnames(all.riskgroups.data[[dataset]][[model]]);
				feature.names <- feature.names[3:length(feature.names)]; # ignore survtime & survstat
				for (feature.name in feature.names) {

					# call appropriate KM plotting function
					eval(
						call(
							KM.plotting.fun,
							riskgroup = all.riskgroups.data[[dataset]][[model]][, feature.name],
							survtime = all.riskgroups.data[[dataset]][[model]][, "survtime"],
							survstat = all.riskgroups.data[[dataset]][[model]][, "survstat"],
							file.name = paste(graphs.dir, "KM_uv__", dataset, "__TRAINING_", all.training.names,  "__model_", model, "__", feature.name,  "__truncated_", truncate.survival, ".png", sep=""),
							main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ")", sep = ""), ""),
							resolution = resolution
							)
						);

					# call sensitivity analysis function
					create.sensitivity.plot(
						riskscore = all.riskscores.data[[dataset]][[model]][, feature.name],
						riskgroup = all.riskgroups.data[[dataset]][[model]][, feature.name],
						survtime = all.riskgroups.data[[dataset]][[model]][, "survtime"],
						survstat = all.riskgroups.data[[dataset]][[model]][, "survstat"],
						survtime.cutoffs = survtime.cutoffs,
						output.directory = output.directory,
						file.stem = paste("sensitivity_analysis_uv__", dataset, "__TRAINING_", all.training.names,  "__model_", model, "__", feature.name, sep=""),
						main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ")", sep = ""), ""),
						resolution = resolution
						);
					}
				}
			}
		}
	}
