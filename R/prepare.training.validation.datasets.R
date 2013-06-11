prepare.training.validation.datasets <- function(data.directory = ".", output.directory = ".", data.types = c("mRNA"), feature.selection.datasets = NULL, datasets = NULL, networks.database = "default", write.normed.datasets = TRUE, subset = NULL) {

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

	# read cancer datasets
	cancer.data <- load.cancer.datasets(
		datasets.to.load = datasets, 
		data.types = data.types, 
		data.directory = data.directory,
		subset = subset
		);

	# lets scale the data
	for (data.type in data.types) {
		for(dataset in names(cancer.data[["all.data"]][[data.type]])) {
			cancer.data[["all.data"]][[data.type]][[dataset]] <- t(scale(t(cancer.data[["all.data"]][[data.type]][[dataset]])));
			}
		}

	# lets write the normalised/scaled datasets to file system
	if (write.normed.datasets == TRUE) {
		for (data.type in data.types) {
			for(dataset in names(cancer.data[["all.data"]][[data.type]])) {

				# write variation data
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

	# lets read coeffients of nodes and edges
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
	p.thresholds <- c(1); # to be extended later to support multiple P values
	for (p.threshold in p.thresholds) { 
		for(dataset in names(cancer.data[["all.survobj"]])) {
			subnet.scores <- NULL; 
			for (subnet in names(all.adjacency.matrices)) {

				# initialise subnet scores datastructure for total number of patients
				init.patient.scores <- rep(0, nrow(cancer.data[["all.survobj"]][[dataset]]));
				subnet.scores[["1"]][[subnet]] <- init.patient.scores;
				subnet.scores[["2"]][[subnet]] <- init.patient.scores;
				subnet.scores[["3"]][[subnet]] <- init.patient.scores;

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
					nodes.valid <- which( coef.nodes.edges[[data.type]][["nodes.coef"]][ nodes, "P" ] <= p.threshold );
					if (length(nodes.valid) > 0) {
						nodes.valid <- nodes[nodes.valid];
						tmp.patient.scores <- coef.nodes.edges[[data.type]][["nodes.coef"]][nodes.valid, "coef"] *
							cancer.data[["all.data"]][[data.type]][[dataset]][nodes.valid, ];

						# handle isolated NA case and Inf/-Inf:
						tmp.patient.scores[is.na(tmp.patient.scores)] <- 0;
						tmp.patient.scores[tmp.patient.scores == Inf | tmp.patient.scores == -Inf] <- 0;

						if (length(nodes.valid) > 1) { tmp.patient.scores <- apply(tmp.patient.scores, 2, sum) };

						subnet.scores[["1"]][[subnet]] <- subnet.scores[["1"]][[subnet]] + tmp.patient.scores;
						subnet.scores[["2"]][[subnet]] <- subnet.scores[["2"]][[subnet]] + tmp.patient.scores;
						}

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
										tmp.patient.scores <- coef.nodes.edges[[data.type]][["edges.coef"]][g1, g2] *
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
					}
				}
			# save patient subnet scores this dataset
			for (model in names(subnet.scores)) {
				x <- t(as.data.frame(subnet.scores[[model]]));
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
					sep="\t"
					);
				}
			}
		}

	#return (cancer.data);
	}
