#' Prepare training and validation datasets
#' 
#' Computes per-patient pathway-derived network impact scores across all input
#' datasets, independently
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files will be saved
#' @param data.types A vector of molecular datatypes to load. Defaults to
#' c('mRNA')
#' @param data.types.ordinal A vector of molecular datatypes to be treated as
#' ordinal. Defaults to c('cna')
#' @param min.ordinal.threshold A named vector specifying minimum percent
#' threshold for each ordinal data type to be used prior to estimating
#' coefficients. Coefficient for features not satisfying minimum threshold will
#' not be estimated, and set to 0. Defaults to cna threshold as 3 percent
#' @param centre.data A character string specifying the centre value to be used for 
#' scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
#' e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
#' as well as for dichotomising data for estimating univariate betas from Cox model.
#' Defaults to 'median'
#' @param p.threshold Cox P value threshold to be applied for selecting features 
#' (e.g. genes) which will contribute to patient risk score estimation. Defaults to 0.5
#' @param feature.selection.datasets A vector containing names of datasets used
#' for feature selection in function \code{derive.network.features()}
#' @param datasets A vector containing names of all the datasets to be later
#' used for training and validation purposes
#' @param truncate.survival A numeric value specifying survival truncation in
#' years. Defaults to 100 years which effectively means no truncation
#' @param networks.database Name of the pathway networks database. Default to
#' NCI PID/Reactome/Biocarta i-e "default"
#' @param write.normed.datasets A toggle to control whether processed mRNA and
#' survival data should be written to file
#' @param subset A list with a Field and Entry component specifying a subset of
#' patients to be selected whose annotation Field matches Entry
#' @return The output files are stored under \code{output.directory}/output/
#' @author Syed Haider
#' @keywords IO
#' @examples
#' 
#' # get data directory 
#' data.directory <- get.program.defaults()[["test.data.dir"]];
#' 
#' # initialise params
#' output.directory <- tempdir();
#' data.types <- c("mRNA");
#' feature.selection.datasets <- c("Breastdata1");
#' training.datasets <- c("Breastdata1");
#' validation.datasets <- c("Breastdata1", "Breastdata2");
#' 
#' # preparing training and validation datasets.
#' # Normalisation & patientwise subnet feature scores
#' prepare.training.validation.datasets(
#'   data.directory = data.directory,
#'   output.directory = output.directory,
#'   data.types =  data.types,
#'   feature.selection.datasets = feature.selection.datasets,
#'   datasets = unique(c(training.datasets, validation.datasets)),
#'   networks.database = "test"
#'   );
#' 
#' @export prepare.training.validation.datasets
prepare.training.validation.datasets <- function(
	data.directory = ".", 
	output.directory = ".", 
	data.types = c("mRNA"), 
	data.types.ordinal = c("cna"), 
	min.ordinal.threshold = c("cna" = 3), 
	centre.data = "median", 
	p.threshold = 0.5, 
	feature.selection.datasets = NULL, 
	datasets = NULL, 
	truncate.survival = 100, 
	networks.database = "default", 
	write.normed.datasets = TRUE, 
	subset = NULL
	) {

	# output directory
	out.dir <- paste(output.directory, "/output/", sep = "");

	# find out where is the data dir bundled with package
	program.data <- get.program.defaults(networks.database = networks.database);

	# program data files and initialise variables
	subnets.file <- program.data[["subnets.file"]];
	all.feature.selection.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.adjacency.matrices <- list();
	subnet.scores <- list();
	cancer.data <- list();
	coef.nodes.edges <- list();
	subnets.selected.features <- list();

	# read cancer datasets
	cancer.data <- load.cancer.datasets(
		truncate.survival = truncate.survival,
		datasets.to.load = datasets, 
		data.types = data.types, 
		data.directory = data.directory,
		subset = subset
		);

	# lets scale the data
	for (data.type in data.types) {
		for(dataset in names(cancer.data[["all.data"]][[data.type]])) {
			if (!(data.type %in% data.types.ordinal)) {
				#cancer.data[["all.data"]][[data.type]][[dataset]] <- 
				#	t(scale(t(cancer.data[["all.data"]][[data.type]][[dataset]])));
				cancer.data[["all.data"]][[data.type]][[dataset]] <- 
					t( centre.scale.dataset(x = t(cancer.data[["all.data"]][[data.type]][[dataset]]), centre.data = centre.data) );
				}
			}
		}

	# lets write the normalised/scaled datasets to file system
	if (write.normed.datasets == TRUE) {
		for (data.type in data.types) {
			for(dataset in names(cancer.data[["all.data"]][[data.type]])) {

				# write molecular data
				write.table(
					x = cancer.data[["all.data"]][[data.type]][[dataset]],
					file = paste(out.dir, "/", dataset, '_', data.type, ".txt", sep = ""),
					sep = "\t",
					col.names = NA
					);

				# write survival data to file (this gets done twice - never mind though)
				y <- cancer.data[["all.survobj"]][[dataset]];
				rownames(y) <- colnames(cancer.data[["all.data"]][[data.type]][[dataset]]);
				cancer.data[["all.survobj"]][[dataset]] <- y;

				write.table(
					x = y,
					file = paste(out.dir, "/", dataset, '_Survival.txt', sep = ''),
					sep = "\t",
					col.names = NA
					);
				}
			}
		}

	# lets read coefficients of nodes and edges
	for (data.type in data.types) {

		# read nodes coefficients & P
		coef.nodes.edges[[data.type]][["nodes.coef"]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/coxph_nodes__", all.feature.selection.names, "__datatype_", data.type, ".txt", sep=""),
				row.names = 1,
				header = TRUE,
				sep = "\t"
				)
			);

		# read nodes' fitted cox models
		con <- gzfile(
			paste(out.dir, "/coxph_nodes__", all.feature.selection.names, "__datatype_", data.type, "_models.gz", sep=""), 
			"rb"
			);
		coef.nodes.edges[[data.type]][["cox.uv"]] <- unserialize(con);
		close(con);

		# read edges coefficients
		coef.nodes.edges[[data.type]][["edges.coef"]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/coxph_edges_coef__", all.feature.selection.names, "__datatype_", data.type, ".txt", sep=""),
				row.names = 1,
				header = TRUE,
				sep = "\t"
				)
			);
		colnames(coef.nodes.edges[[data.type]][["edges.coef"]]) <- gsub("^X", "", colnames(coef.nodes.edges[[data.type]][["edges.coef"]]), perl = TRUE);

		# read edges P
		coef.nodes.edges[[data.type]][["edges.P"]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/coxph_edges_P__", all.feature.selection.names, "__datatype_", data.type, ".txt", sep=""),
				row.names = 1,
				header = TRUE,
				sep = "\t"
				)
			);
		colnames(coef.nodes.edges[[data.type]][["edges.P"]]) <- gsub("^X", "", colnames(coef.nodes.edges[[data.type]][["edges.P"]]), perl = TRUE);
		}

	# lets read input subnetworks
	all.adjacency.matrices <- get.adjacency.matrix(subnets.file);

	# add "_at" to entrez ids
	for (subnet in names(all.adjacency.matrices)) {
		rownames(all.adjacency.matrices[[subnet]]) <- paste(rownames(all.adjacency.matrices[[subnet]]), "_at", sep="");
		colnames(all.adjacency.matrices[[subnet]]) <- paste(colnames(all.adjacency.matrices[[subnet]]), "_at", sep="");
		}

	# lets compute patient-wise subnetwork scores
	# to be extended later to support multiple P values,
	for(dataset in names(cancer.data[["all.survobj"]])) {
		subnet.scores <- NULL; 
		for (subnet in names(all.adjacency.matrices)) {

			# initialise subnet scores datastructure for total number of patients
			init.patient.scores <- rep(0, nrow(cancer.data[["all.survobj"]][[dataset]]));

			subnet.scores[["1"]][[subnet]] <- init.patient.scores;
			subnet.scores[["2"]][[subnet]] <- init.patient.scores;
			subnet.scores[["3"]][[subnet]] <- init.patient.scores;

			# track nodes and interactions that are seen already by one data-type
			nodes.seen <- edges.seen <- NULL;

			for(data.type in data.types) {
				adjacency.matrix <- all.adjacency.matrices[[subnet]];
				nodes <- rownames(adjacency.matrix);
				nodes <- intersect(
					nodes,
					intersect(
						rownames(coef.nodes.edges[[data.type]][["nodes.coef"]]),
						rownames(cancer.data[["all.data"]][[data.type]][[dataset]])
						)
					);
				adjacency.matrix <- adjacency.matrix[nodes, nodes];

				# Nodes score
				# NOTE: for ordinal data types, feature selection is conducted on the 
				# beta/coef for smallest P value group/level. However, predict() calculates
				# the patientwise scores correctly given the group/level
				nodes.valid <- which( coef.nodes.edges[[data.type]][["nodes.coef"]][ nodes, "P" ] <= p.threshold );
				nodes.valid <- nodes[nodes.valid];

				# store features to save to file system 
				subnets.selected.features[[data.type]][[subnet]] <- vector();
				subnets.selected.features[[data.type]][[subnet]] <- nodes.valid;

				# [inactivated] track nodes and interactions that are seen already by one data-type
				# and therefore be ignored by the second one if already in the model.
				# Respects the data.types order eg. mRNA, cna would give mRNA priority
				# nodes.valid <- setdiff(nodes.valid, nodes.seen);
				# nodes.seen <- union(nodes.seen, nodes.valid);

				if (length(nodes.valid) > 0) {

					# initialise variables
					tmp.patient.scores <- NULL;

					# process ordinal variables differently
					if (data.type %in% data.types.ordinal) {

						# find minimum ordinal threshold for this datatype
						ordinal.threshold <- length(init.patient.scores) * 
							min.ordinal.threshold[data.type]/100;

						# go through feature by feature and calculate per-patient score
						tmp.patient.scores <- sapply(
							nodes.valid, 
							FUN = function(node.valid) {

								coef.vector <- coef.nodes.edges[[data.type]][["cox.uv"]][[node.valid]][, "coef"];
								names(coef.vector) <- as.character(rownames(coef.nodes.edges[[data.type]][["cox.uv"]][[node.valid]]));

								# set coef of levels with p > threshold to 0, this is often the case resulting in large betas
								p.vector <- coef.nodes.edges[[data.type]][["cox.uv"]][[node.valid]][, "p"];
								which.p.invalid <- which(p.vector > p.threshold);
								if (length(which.p.invalid) > 0) {
									coef.vector[which.p.invalid] <- 0;
									}

								# ordinal data without a particular level can have beta = NA as well. set to 0
								coef.vector[is.na(coef.vector)] <- 0;

								# set coef vector to 0 if less than X percent had aberrations as this may
								# result in large betas
								if (max(
									coef.nodes.edges[[data.type]][["cox.uv"]][[node.valid]][, "n"], 
									na.rm = TRUE
									) < ordinal.threshold) {

										coef.vector[names(coef.vector)] <- 0;

										# remove from selected genes
										subnets.selected.features[[data.type]][[subnet]] <<- subnets.selected.features[[data.type]][[subnet]][
											-which(subnets.selected.features[[data.type]][[subnet]] == node.valid)
											]; 
									}

								# m <- model.matrix(~factor(
								#	t(cancer.data[["all.data"]][[data.type]][[dataset]][node.valid, ]), 
								#	levels = c("0", "-1", "1")
								#	), data = lung)[,][, -1, drop = FALSE];
								# m %*% as.matrix(coef.vector)
								# below is same but a generic solution to whats above where we 
								# have to specify the levels correctly.
								# However, predict() in R is slightly different as it also 
								# subtract means which is not necessary given that we dont do that for 
								# continuous variables either and make linear predictions instead
								unlist(as.vector(
									coef.vector[as.character(cancer.data[["all.data"]][[data.type]][[dataset]][node.valid, ])]
									));
								}
							);

						# transpose or vectorise
						if (length(nodes.valid) > 1) {
							tmp.patient.scores <- t(tmp.patient.scores);
							}
						else {
							tmp.patient.scores <- tmp.patient.scores[, 1];
							}
						}
					else {
						tmp.patient.scores <- coef.nodes.edges[[data.type]][["nodes.coef"]][nodes.valid, "coef"] *
							cancer.data[["all.data"]][[data.type]][[dataset]][nodes.valid, ];
						}

					# handle isolated NA case and Inf/-Inf:
					tmp.patient.scores[is.na(tmp.patient.scores)] <- 0;
					tmp.patient.scores[tmp.patient.scores == Inf | tmp.patient.scores == -Inf] <- 0;

					if (length(nodes.valid) > 1) { tmp.patient.scores <- apply(tmp.patient.scores, 2, sum) };

					subnet.scores[["1"]][[subnet]] <- subnet.scores[["1"]][[subnet]] + tmp.patient.scores;
					subnet.scores[["2"]][[subnet]] <- subnet.scores[["2"]][[subnet]] + tmp.patient.scores;
					}

				# recover memory
				gc();

				# Edges score
				if (length(nodes) > 1) { # ignore self interaction matrix
					for (i in 2:nrow(adjacency.matrix)) {
						j.max <- i - 1;
						edge.cols <- which(adjacency.matrix[i, 1:j.max] == 1);
						if (length(edge.cols) > 0) {
							g1 <- nodes[i];
							for (nodes.i in 1:length(edge.cols)) {
								g2 <- nodes[edge.cols[nodes.i]];
								if (!is.na(coef.nodes.edges[[data.type]][["edges.P"]][g1, g2]) 
									& coef.nodes.edges[[data.type]][["edges.P"]][g1, g2] <= p.threshold) {
									tmp.patient.scores <- abs(coef.nodes.edges[[data.type]][["edges.coef"]][g1, g2]) *
										cancer.data[["all.data"]][[data.type]][[dataset]][g1, ] *
										cancer.data[["all.data"]][[data.type]][[dataset]][g2, ];

									# handle isolated NA case and Inf/-Inf:
									tmp.patient.scores[is.na(tmp.patient.scores)] <- 0;
									tmp.patient.scores[tmp.patient.scores == Inf | tmp.patient.scores == -Inf] <- 0;

									subnet.scores[["1"]][[subnet]] <- subnet.scores[["1"]][[subnet]] + tmp.patient.scores;
									subnet.scores[["3"]][[subnet]] <- subnet.scores[["3"]][[subnet]] + tmp.patient.scores;
									}
								}
							}
						}
					}

				# collapse the list of selected genes by comma
				subnets.selected.features[[data.type]][[subnet]] <- paste(
					subnets.selected.features[[data.type]][[subnet]], collapse = ", "
					);
				}
			}

		# save patient subnet scores for this dataset
		for (model in names(subnet.scores)) {
			x <- do.call(rbind, subnet.scores[[model]]);
			x <- rbind(
				"survtime" = cancer.data[["all.survobj"]][[dataset]][, "time"],
				"survstat" = cancer.data[["all.survobj"]][[dataset]][, "status"],
				x
				);
			write.table(
				x = x,
				file = paste(out.dir, "/patientwise_subnets_score__", dataset, "__TRAINING_", all.feature.selection.names, "__model_", model, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);
			}

		# save selected features per subnet for this dataset
		for (data.type in names(subnets.selected.features)) {
			write.table(
				x = as.matrix(subnets.selected.features[[data.type]]),
				file = paste(out.dir, "/subnets_selected_features__", dataset, "__TRAINING_", all.feature.selection.names, "__datatype_", data.type, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);
			}

		} # end of dataset loop

	# recover memory
	gc.tmp <- gc();

	# return ();
	}
