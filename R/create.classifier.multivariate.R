#' Trains and tests a multivariate survival model
#' 
#' Trains a model on training datasets. Predicts the risk score for all the
#' training & datasets, independently. This function also predicts the risk
#' score for combined training datasets cohort and validation datasets cohort.
#' The risk score estimation is done by multivariate models fit by
#' \code{fit.survivalmodel}. The function also predicts risk scores for each of
#' the \code{top.n.features} independently.
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{feature.selection.datasets}, \code{training.datasets},
#' \code{validation.datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files will be saved
#' @param feature.selection.datasets A vector containing names of datasets used
#' for feature selection in function \code{derive.network.features()}
#' @param feature.selection.p.threshold One of the P values that were used for
#' feature selection in function \code{derive.network.features()}. This
#' function does not support vector of P values as used in
#' \code{derive.network.features()} for performance reasons
#' @param training.datasets A vector containing names of training datasets
#' @param validation.datasets A vector containing names of validation datasets
#' @param top.n.features A numeric value specifying how many top ranked
#' features will be used for univariate survival modelling
#' @param models A character vector specifying which of the models ('1' = N+E,
#' '2' = N, '3' = E) to run
#' @param learning.algorithms A character vector specifying which learning
#' algorithm to be used for model fitting and feature selection. Defaults to
#' c('backward', 'forward'). Available options are: c('backward', 'forward',
#' 'glm', 'randomforest')
#' @param alpha.glm A numeric vector specifying elastic-net mixing parameter
#' alpha, with range alpha raning from [0,1]. 1 for LASSO (default) and 0 for
#' ridge. For multiple values of alpha, most optimal value is selected through
#' cross validation on training set
#' @param k.fold.glm A numeric value specifying k-fold cross validation if glm
#' was chosen in \code{learning.algorithms}
#' @param seed.value A numeric value specifying seed for glm k-fold cross or random forest
#' validation if glm was chosen in \code{learning.algorithms}
#' @param cores.glm An integer value specifying number of cores to be used for
#' glm if it was chosen in \code{learning.algorithms}
#' @param rf.ntree An integer value specifying the number of trees in random forest. 
#' Defaults to 1000. This should be tuned after starting with a large forest such as 1000 
#' in the initial run and assessing the results in output\/OOB_error__TRAINING_* to see 
#' where the OOB error rate stablises, and then rerunning with the stablised rf.ntree 
#' parameter
#' @param rf.mtry An integer value specifying the number of variables randomly selected 
#' for splitting a node. Defaults to sqrt(features), which is the same as in the 
#' underlying R package random survival forest \code{randomForestSRC::rfsrc} 
#' @param rf.nodesize An integer value specifying number of unique cases in a terminal 
#' node. Defaults to 15, which is the same as in the underlying R package random survival 
#' forest \code{randomForestSRC::rfsrc} 
#' @param rf.samptype An character string specifying name of sampling. Defaults to 
#' sampling without replacement 'swor'. Available options are: c('swor', 'swr')
#' @param rf.sampsize A function specifying sampling size when \code{rf.samptype} is set 
#' to sampling without replacement ('swor'). Defaults to 66\%: \code{function(x){x * .66}}
#' @param ... other params to be passed on to the random forest call to the underlying 
#' R package random survival forest \code{randomForestSRC::rfsrc} 
#' 
#' @return The output files are stored under \code{output.directory}/output/
#' @author Syed Haider & Vincent Stimper
#' @keywords survival
#' @examples
#' 
#' # see package's main documentation
#' 
#' @import randomForestSRC
#' @export create.classifier.multivariate
create.classifier.multivariate <- function(
	data.directory = ".", 
	output.directory = ".", 
	feature.selection.datasets = NULL, 
	feature.selection.p.threshold = 0.05, 
	training.datasets = NULL, 
	validation.datasets = NULL, 
	top.n.features = 25, 
	models = c("1", "2", "3"), 
	learning.algorithms = c("backward", "forward"), 
	alpha.glm = c(1), 
	k.fold.glm = 10, 
	seed.value = 51214, 
	cores.glm = 1,
	rf.ntree = 1000,
	rf.mtry = NULL,
	rf.nodesize = 15,
	rf.samptype = "swor",
	rf.sampsize = function(x){x * .66},
	...
	) {

	# sanity checks
	if (top.n.features < 2) {
		cat("\ntop.n.features < 2, hence no need to fit a multivariate model. Univariate model results should suffice");
		return ();
		}

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# program data files and initialise variables
	all.feature.selection.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.training.names <- paste(sort(training.datasets), collapse="_");
	top.subnets <- list();
	all.subnet.scores <- list();
	all.training.data <- NULL;
	all.validation.data <- NULL;
	all.training.median <- NULL;
	stepAIC.results <- NULL;
	coxph.header <- c("HR", "95l", "95u", "P", "n");

	# lets read top subnets for the give feature.selection.datasets
	for (model in models) {

		# read in all subnets scores
		top.subnets[[model]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/top_subnets_score__TRAINING_", all.feature.selection.names, "__model_", model, "__PV_", feature.selection.p.threshold, ".txt", sep = ""),
				header = TRUE,
				row.names = 1,
				sep = "\t"
				)
			);

		# read training and validation datasets
		for (dataset in c(training.datasets, validation.datasets)) {

			all.subnet.scores[[dataset]][[model]] <- as.matrix(
				read.table(
					file = paste(out.dir, "patientwise_subnets_score__", dataset, "__TRAINING_", all.feature.selection.names, "__model_", model, ".txt", sep = ""),
					header = TRUE,
					row.names = 1,
					sep = "\t"
					)
				);

			# join the training & validation datasets together to make two large independent training & validation datasets
			if (dataset %in% training.datasets) {
				all.training.data[[model]] <- cbind(all.training.data[[model]], all.subnet.scores[[dataset]][[model]]);
				}
			if (dataset %in% validation.datasets) {
				all.validation.data[[model]] <- cbind(all.validation.data[[model]], all.subnet.scores[[dataset]][[model]]);
				}
			}
		}
	
	# create and register cluster in case of parallel computing for glmnet
	if ("glm" %in% learning.algorithms) {
		if (cores.glm > 1) {
			computing.cluster <- makeCluster(cores.glm);
			registerDoParallel(computing.cluster);
			} else {
			registerDoSEQ();
			}
		}

	# run learning algorithm for each model, for feature selection followed by patient risk score prediction
	for (model in models) {

		# in case requested number of features are greater than available features
		if (top.n.features > nrow(top.subnets[[model]])) {
			stop("\ntop.n.features larger than available features i-e ", nrow(top.subnets[[model]]));
			}

		explanatory.variables <- rownames(top.subnets[[model]])[1:top.n.features];
		training.data <- data.frame(
			t(all.training.data[[model]][c("survtime", "survstat", explanatory.variables), ])
			);
		validation.data <- data.frame(
			t(all.validation.data[[model]][c("survtime", "survstat", explanatory.variables), ])
			);

		for (learning.algorithm in learning.algorithms) {

			cat("\n*****************  starting fit\tModel ", model, learning.algorithm, "  *****************\n");
			model.formula <- NULL;
			model.fit <- NULL;
			selected.features <- vector();

			# Forward/backward refinement algorithms
			if (learning.algorithm %in% c("backward", "forward")) {
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

				print(model.formula[[learning.algorithm]]);
				tryCatch(
					expr = {
						model.fit <- stepAIC(
							coxph(model.formula[[learning.algorithm]], data = training.data),
							direction = learning.algorithm,
							trace = TRUE,
							scope = list(
								lower = coxph(model.formula[["forward"]], data = training.data),
								upper = coxph(model.formula[["backward"]], data = training.data) 
								)
							)
						},
					error = function(ex) {
						cat("\nModel failed to converge (a known coxph issue) after stepAIC using MDS model: ", model);
						}
					);

				# prepare model fit's summary
				model.summary <- summary(model.fit);

				# break if necessary
				# model did not select any variable
				if (!("conf.int" %in% names(model.summary)) ||
					"conf.int" %in% names(model.summary) && is.null(model.summary$conf.int)) {
					cat("\nNull Model, no variables selected after stepAIC using MDS model: ", model);
					next;
					}
				# fail to converge case
				if ("Class" %in% names(model.summary) && "Mode" %in% names(model.summary) 
					&& model.summary[["Class"]] == "NULL" && model.summary[["Mode"]] == "NULL") {
					next;
					}

				# store the results of each model (summary) to file
				res.matrix <- matrix(
					data = NA,
					ncol = 5,
					nrow = nrow(model.summary$conf.int),
					dimnames = list(
						rownames(model.summary$conf.int),
						c("HR", "95l", "95u", "p-val", "beta")
						)
					);
				selected.features <- rownames(res.matrix);

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
					file = paste(out.dir, "/beta__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				cat("\nstoring AIC TRACE, the last one being the AIC of the final model");
				write.table(
					as.matrix(model.fit$anova[1:6]),
					file = paste(out.dir, "/AIC__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			# GLM net refinement algorithm
			if (learning.algorithm == "glm") {

				which.surv.cols <- which(colnames(training.data) %in% c("survtime", "survstat"));
				x.T <- as.matrix(training.data[, -which.surv.cols]);
				y.T <- cbind("time" = training.data[, "survtime"], "status" = training.data[, "survstat"]);

				# initialise glm params
				set.seed(seed.value);
				foldid <- sample(rep(seq(k.fold.glm), length.out = nrow(x.T)));
				cvm <- rep(NA, length(alpha.glm));
				model.fit.tmp <- vector("list", length(alpha.glm));

				# traverse through all alpha and perform cross validation
				if (length(alpha.glm) > k.fold.glm) {
					model.fit.tmp <- foreach(i = seq_along(alpha.glm), .packages = c("survival", "glmnet")) %dopar% {
						cv.glmnet(
							x = x.T, 
							y = y.T, 
							alpha = alpha.glm[i], 
							foldid = foldid, 
							standardize = FALSE, 
							family = "cox"
							)
						}
					} else {
					for (i in seq_along(alpha.glm)) {
						model.fit.tmp[[i]] <- cv.glmnet(
							x = x.T, 
							y = y.T, 
							alpha = alpha.glm[i], 
							foldid = foldid, 
							standardize = FALSE, 
							family = "cox",
							parallel = TRUE
							);
						}
					}

				# extract cvm for model with best (min) lambda
				cvm <- sapply(model.fit.tmp, function(mft) mft$cvm[which(mft$lambda == mft$lambda.min)]);

				# extract the optimal model
				model.fit <- model.fit.tmp[[which.min(cvm)]];
				remove(model.fit.tmp);

				# select optimal lambda
				lambda.to.use <- model.fit$lambda.min;

				# features selected by the final model
				selected.features <- names(which(model.fit$glmnet.fit$beta[, 
					which(model.fit$lambda == lambda.to.use)
					] != 0));

				if (length(selected.features) == 0) {
					cat("\nNull Model, no variables selected after GLM using MDS model: ", model);
					next;
					}

				glm.coef <- as.matrix(coef(model.fit, s = "lambda.min"));
				colnames(glm.coef) <- "beta";

				cat("\nstoring GLM fit for: ", model);
				write.table(
					glm.coef[selected.features, , drop = FALSE],
					file = paste(out.dir, "/beta__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			# Random survival forest implementation
			if (learning.algorithm == "randomforest") {

				# create random forest
				model.fit <- rfsrc(
					Surv(survtime, survstat) ~ ., 
					data = training.data, 
					importance = TRUE,
					ntree = rf.ntree,
					mtry = ifelse(is.null(rf.mtry), round(sqrt(ncol(training.data)-2)), rf.mtry),
					nodesize = rf.nodesize,
					samptype = rf.samptype,
					sampsize = rf.sampsize,
					block.size = 1, # calculate error rate on every nth tree. default to NULL with error rate of last tree only,
					var.used = "by.tree",
					statistics = TRUE,
					membership = TRUE,
					seed = seed.value
					);

				print(model.fit);

				# store random forest object for subsequent post-processing of trees
				save(
					model.fit, 
					file = paste(out.dir, "/model_object__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".RData", sep = "")
					);

				# store random forest vimp
				cat("\nstoring Random Forest variables with their importance");
				write.table(
					cbind("vimp" = sort(model.fit$importance, decreasing = TRUE)),
					file = paste(out.dir, "/vimp__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				# store random forest OOB error rate
				cat("\nstoring Random Forest OOB error rate");
				write.table(
					cbind("Number of trees" = 1:length(model.fit$err.rate), "OOB error rate" = model.fit$err.rate),
					file = paste(out.dir, "/OOB_error__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					row.names = FALSE,
					sep = "\t"
					);

				}

			# estimate risk scores over TRAINING datasets
			all.training.risk.scores <- NULL;
			all.training.groups <- NULL;
			all.training.survtime <- NULL;
			all.training.survstat <- NULL;

			for (dataset in training.datasets) {

				# prepare data for prediction
				new.data <- data.frame(t(all.subnet.scores[[dataset]][[model]]));
				new.data <- new.data[, c("survtime", "survstat", explanatory.variables)];

				# predict per-patient risk score
				if (learning.algorithm %in% c("backward", "forward")) {
					risk.scores <- predict(model.fit, newdata = new.data, type = "risk");
					}
				if (learning.algorithm == "glm") {
					risk.scores <- predict(
						model.fit, 
						newx = as.matrix(new.data[, -which(colnames(new.data) %in% c("survtime", "survstat"))]), 
						type = "response", 
						s = lambda.to.use
						)[, 1];
					}
				if (learning.algorithm == "randomforest") {
					risk.scores <- predict(
						model.fit, 
						newdata = new.data[, -which(colnames(new.data) %in% c("survtime", "survstat"))] 
						)$predicted;
					}

				risk.groups <- SIMMS::dichotomize.dataset(risk.scores);
				names(risk.scores) <- rownames(new.data);
				names(risk.groups) <- rownames(new.data);
				survtime <- new.data$survtime;
				survstat <- new.data$survstat;

				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				names(coxph.res) <- coxph.header;

				# save the results to file system for later visualizations (either BL or ordinary graphics)
				write.table(
					x = cbind("riskgroup" = risk.groups, "survtime" = survtime, "survstat" = survstat),
					file = paste(out.dir, "riskgroups__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = cbind("riskscore" = risk.scores, "survtime" = survtime, "survstat" = survstat),
					file = paste(out.dir, "riskscores__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = t(coxph.res),
					file = paste(out.dir, "coxph__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				if (!is.na(coxph.res[4])) {
					survdiff.res <- get.chisq.stats(risk.groups, Surv(survtime, survstat));
					write.table(
						x = t(survdiff.res),
						file = paste(out.dir, "survdiff__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}

				all.training.risk.scores <- c(all.training.risk.scores, risk.scores);
				all.training.survtime <- c(all.training.survtime, survtime);
				all.training.survstat <- c(all.training.survstat, survstat);
				}

			# combined training set
			all.training.groups <- SIMMS::dichotomize.dataset(all.training.risk.scores);
			names(all.training.groups) <- names(all.training.risk.scores);
			coxph.res <- SIMMS::fit.coxmodel(
				all.training.groups, Surv(all.training.survtime, all.training.survstat)
				)$cox.stats;
			names(coxph.res) <- coxph.header;
	
			write.table(
				x = cbind("riskgroup" = all.training.groups, "survtime" = all.training.survtime, "survstat" = all.training.survstat),
				file = paste(out.dir, "riskgroups__", "all_training", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind("riskscore" = all.training.risk.scores, "survtime" = all.training.survtime, "survstat" = all.training.survstat),
				file = paste(out.dir, "riskscores__", "all_training", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = t(coxph.res),
				file = paste(out.dir, "coxph__", "all_training", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			if (!is.na(coxph.res[4])) {
				survdiff.res <- get.chisq.stats(
					all.training.groups, Surv(all.training.survtime, all.training.survstat)
					);
				write.table(
					x = t(survdiff.res),
					file = paste(out.dir, "survdiff__", "all_training", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			all.training.median <- median(all.training.risk.scores);

			# store all_training set median score
			write(
				x = all.training.median,
				file =  paste(out.dir, "riskscore_median__", "all_training", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep="")
				);
			
			# estimate risk scores over VALIDATION datasets
			all.validation.risk.scores <- NULL;
			all.validation.groups <- NULL;
			all.validation.survtime <- NULL;
			all.validation.survstat <- NULL;

			for (dataset in validation.datasets) {

				# prepare data for prediction
				new.data <- data.frame(t(all.subnet.scores[[dataset]][[model]]));
				new.data <- new.data[, c("survtime", "survstat", explanatory.variables)];

				# predict per-patient risk score
				if (learning.algorithm %in% c("backward", "forward")) {
					risk.scores <- predict(model.fit, newdata = new.data, type = "risk");
					}
				if (learning.algorithm %in% c("glm")) {
					lambda.to.use <- model.fit$lambda.min;
					risk.scores <- predict(
						model.fit, 
						newx = as.matrix(new.data[, -which(colnames(new.data) %in% c("survtime", "survstat"))]), 
						type = "response", 
						s = lambda.to.use
						)[, 1];
					}
				if (learning.algorithm == "randomforest") {
					risk.scores <- predict(
						model.fit, 
						newdata = new.data[, -which(colnames(new.data) %in% c("survtime", "survstat"))]
						)$predicted;
					}

				risk.groups <- SIMMS::dichotomize.dataset(risk.scores, split.at = all.training.median);
				names(risk.scores) <- rownames(new.data); 
				names(risk.groups) <- rownames(new.data); 
				survtime <- new.data$survtime;
				survstat <- new.data$survstat;

				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				names(coxph.res) <- coxph.header;

				# save the results to file system for later visualizations (either BL or ordinary graphics)
				write.table(
					x = cbind("riskgroup" = risk.groups, "survtime" = survtime, "survstat" = survstat),
					file = paste(out.dir, "riskgroups__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = cbind("riskscore" = risk.scores, "survtime" = survtime, "survstat" = survstat),
					file = paste(out.dir, "riskscores__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = t(coxph.res),
					file = paste(out.dir, "coxph__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				
				if (!is.na(coxph.res[4])) {
					survdiff.res <- get.chisq.stats(risk.groups, Surv(survtime, survstat));
					write.table(
						x = t(survdiff.res),
						file = paste(out.dir, "survdiff__", dataset, "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}
				
				all.validation.risk.scores <- c(all.validation.risk.scores, risk.scores);
				all.validation.survtime <- c(all.validation.survtime, survtime);
				all.validation.survstat <- c(all.validation.survstat, survstat);
				}

			# combined validation set
			all.validation.groups <- SIMMS::dichotomize.dataset(all.validation.risk.scores, split.at = all.training.median);
			names(all.validation.groups) <- names(all.validation.risk.scores);
			coxph.res <- SIMMS::fit.coxmodel(
				all.validation.groups, Surv(all.validation.survtime, all.validation.survstat)
				)$cox.stats;
			names(coxph.res) <- coxph.header;

			write.table(
				x = cbind("riskgroup" = all.validation.groups, "survtime" = all.validation.survtime, "survstat" = all.validation.survstat),
				file = paste(out.dir, "riskgroups__", "all_validation", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind("riskscore" = all.validation.risk.scores, "survtime" = all.validation.survtime, "survstat" = all.validation.survstat),
				file = paste(out.dir, "riskscores__", "all_validation", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = t(coxph.res),
				file = paste(out.dir, "coxph__", "all_validation", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			if (!is.na(coxph.res[4])) {
				survdiff.res <- get.chisq.stats(
					all.validation.groups, Surv(all.validation.survtime, all.validation.survstat)
					);
				write.table(
					x = t(survdiff.res),
					file = paste(out.dir, "survdiff__", "all_validation", "__TRAINING_", all.training.names, "__", learning.algorithm, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}
			}
	  }
	  
  	# stop cluster
	if ("glm" %in% learning.algorithms) {
		if (cores.glm > 1) {
			registerDoSEQ();
			stopCluster(computing.cluster);
			remove(computing.cluster);
			}
		}
	}
