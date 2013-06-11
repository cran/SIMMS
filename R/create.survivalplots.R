create.survivalplots <- function(data.directory = ".", output.directory = ".", training.datasets = NULL, validation.datasets = NULL, top.n.features = 25, truncate.survival = 100, survtime.cutoffs = c(seq(5,10,1)), main.title = FALSE, KM.plotting.fun = "create.KM.plot", resolution = 100) {

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");
	logs.dir <- paste(output.directory, "/logs/", sep = "");

	# program data files and initialise variables
	all.training.names <- paste(sort(training.datasets), collapse="_");
	all.validation.names <- paste(sort(validation.datasets), collapse="_");
	models <- c("1", "2", "3");
	model.names <- c("N+E", "N", "E");
	all.riskgroups.data <- list();
	all.riskscores.data <- list();

	# lets read in the training & validation dataset riskgroups (stepAIC selected multivariate model)
	for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
		for (model in models) {
			for (direction in c("forward", "backward")) {

				all.riskgroups.data[[dataset]][[model]][[direction]] <- as.matrix(
					read.table(
						file = paste(out.dir, "/riskgroups__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						header = TRUE,
						row.names = 1,
						sep = "\t"
						)
					);

				all.riskscores.data[[dataset]][[model]][[direction]] <- as.matrix(
					read.table(
						file = paste(out.dir, "/riskscores__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						header = TRUE,
						row.names = 1,
						sep = "\t"
						)
					);
				}
			}
		}

	# lets do KM plotting (stepAIC selected multivariate model)
	for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
		for (model in models) {
			for (direction in c("forward", "backward")) {

				# call appropriate KM plotting function, user defined OR this package's default
				eval(
					call(
						KM.plotting.fun,
						riskgroup = all.riskgroups.data[[dataset]][[model]][[direction]][, "riskgroup"],
						survtime = all.riskgroups.data[[dataset]][[model]][[direction]][, "survtime"],
						survstat = all.riskgroups.data[[dataset]][[model]][[direction]][, "survstat"],
						truncate.survival = truncate.survival,
						file.name = paste(graphs.dir, "KM__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, "__truncated_", truncate.survival, ".tiff", sep=""),
						main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ", ", direction, ")", sep = ""), ""),
						resolution = resolution
						)
					);

				# call sensitivity analysis function
				create.sensitivity.plot(
					riskscore = all.riskscores.data[[dataset]][[model]][[direction]][, "riskscore"],
					riskgroup = all.riskgroups.data[[dataset]][[model]][[direction]][, "riskgroup"],
					survtime = all.riskgroups.data[[dataset]][[model]][[direction]][, "survtime"],
					survstat = all.riskgroups.data[[dataset]][[model]][[direction]][, "survstat"],
					survtime.cutoffs = survtime.cutoffs,
					output.directory = output.directory,
					file.stem = paste("sensitivity_analysis__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, sep=""),
					main.title = ifelse(main.title, paste(dataset, "(", model.names[as.numeric(model)], ", ", direction, ")", sep = ""), ""),
					resolution = resolution
					);
				}
			}
		}

	# lets read in the training & validation dataset riskgroups (univariate grouping)
	all.riskgroups.data <- list();
	all.riskscores.data <- list();

	for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
		for (model in models) {

			all.riskgroups.data[[dataset]][[model]] <- as.matrix(
				read.table(
					file = paste(out.dir, "riskgroups_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					header = TRUE,
					row.names = 1,
					sep = "\t"
					)
				);

			all.riskscores.data[[dataset]][[model]] <- as.matrix(
				read.table(
					file = paste(out.dir, "riskscores_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					header = TRUE,
					row.names = 1,
					sep = "\t"
					)
				);
			}
		}

	# let's do KM plotting (univariate grouping)
	for (dataset in c(training.datasets, "all_training", validation.datasets, "all_validation")) {
		for (model in models) {
			feature.names <- colnames(all.riskgroups.data[[dataset]][[model]]);
			feature.names <- feature.names[3:length(feature.names)]; # ignore survtime & survstat
			for (feature.name in feature.names) {

				# call appropriate KM plotting function
				eval(
					call(
						KM.plotting.fun,
						riskgroup = all.riskgroups.data[[dataset]][[model]][, feature.name],
						survtime = all.riskgroups.data[[dataset]][[model]][, "survtime"],
						survstat = all.riskgroups.data[[dataset]][[model]][, "survstat"],
						truncate.survival,
						file.name = paste(graphs.dir, "KM_uv__", dataset, "__TRAINING_", all.training.names,  "__model_", model, "__", feature.name,  "__truncated_", truncate.survival, ".tiff", sep=""),
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
