#' Derive univariate features from pathway-derived networks
#' 
#' This function fits Cox model to features as well as interaction between
#' features. The coefficients of features are subsequently used to compute
#' impact score of each of the pathway-derived networks.
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{feature.selection.datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files will be saved
#' @param data.types A vector of molecular datatypes to load. Defaults to
#' c('mRNA')
#' @param data.types.ordinal A vector of molecular datatypes to be treated as
#' ordinal. Defaults to c('cna')
#' @param feature.selection.fun Name of the function to be used to estimate
#' network coefficients. Defaults to 'calculate.network.coefficients'
#' @param feature.selection.datasets A vector containing names of training
#' datasets to be used to compute cox statistics
#' @param feature.selection.p.thresholds A vector containing P values to be
#' used as threshold for including features into overall impact score of a
#' network
#' @param truncate.survival A numeric value specifying survival truncation in
#' years. Defaults to 100 years which effectively means no truncation
#' @param networks.database Name of the pathway networks database. Default to
#' NCI PID/Reactome/Biocarta i-e "default"
#' @param subset A list with a Field and Entry component specifying a subset of
#' patients to be selected from each dataset whose annotation Field matches
#' Entry
#' @param ... other params to be passed on to user-defined method for
#' estimating coefficients of network features
#' @return The output files are stored under \code{data.directory}/output/
#' @author Syed Haider
#' @keywords FeatureSelection
#' @examples
#' 
#' options("warn" = -1);
#' 
#' # get data directory 
#' data.directory <- get.program.defaults(networks.database = "test")[["test.data.dir"]];
#' 
#' # initialise params
#' output.directory <- ".";
#' data.types <- c("mRNA");
#' feature.selection.datasets <- c("Breastdata1");
#' feature.selection.p.thresholds <- c(0.05);
#' 
#' # estimate network coefficients for all the subnet features
#' derive.network.features(
#'   data.directory = data.directory,
#'   output.directory = output.directory,
#'   data.types = data.types,
#'   feature.selection.fun = "calculate.network.coefficients",
#'   feature.selection.datasets = feature.selection.datasets,
#'   feature.selection.p.thresholds = feature.selection.p.thresholds,
#'   networks.database = "test"
#'   );
#' 
#' @export derive.network.features
derive.network.features <- function(data.directory = ".", output.directory = ".", data.types = c("mRNA"), data.types.ordinal = c("cna"), feature.selection.fun = "calculate.network.coefficients", feature.selection.datasets = NULL, feature.selection.p.thresholds = c(0.05), truncate.survival = 100, networks.database = "default", subset = NULL, ...) {

	# verify that we got appropriate input data
	to.abort <- FALSE;

	# stop processing if we have bad data
	if (to.abort) {
		stop("Inputs failed sanity-checking");
		}

	# find out where is the data dir bundled with package containing networks database
	program.data <- get.program.defaults(networks.database = networks.database);

	# program data files and initialise variables
	subnets.file <- program.data[["subnets.file"]];
	subnets.file.flattened <- program.data[["subnets.file.flattened"]];
	all.training.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.adjacency.matrices <- list();
	subnet.scores <- list();

	# create the output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	dir.create(out.dir, recursive = TRUE);
	dir.create(graphs.dir, recursive = TRUE);

	# create a function dynamically
	dynamic.function <- get(feature.selection.fun, mode = "function");
	coef.nodes.edges <- dynamic.function(
		data.directory = data.directory,
		output.directory = out.dir,
		training.datasets = feature.selection.datasets,
		data.types = data.types,
		data.types.ordinal = data.types.ordinal,
		subnets.file.flattened = subnets.file.flattened,
		truncate.survival = truncate.survival,
		subset = subset,
		...
		);

	# write coefficients to the filesystem
	for(data.type in data.types) {

		write.table(
			x = coef.nodes.edges[[data.type]][["nodes.coef"]],
			file = paste(out.dir, "/coxph_nodes__", all.training.names, "__datatype_", data.type, ".txt", sep=""),
			row.names = TRUE,
			col.names = NA,
			sep = "\t"
			);

		write.table(
			x = coef.nodes.edges[[data.type]][["edges.coef"]],
			file = paste(out.dir, "/coxph_edges_coef__", all.training.names, "__datatype_", data.type, ".txt", sep=""),
			row.names = TRUE,
			col.names = NA,
			sep = "\t"
			);

		write.table(
			x = coef.nodes.edges[[data.type]][["edges.P"]],
			file = paste(out.dir, "/coxph_edges_P__", all.training.names, "__datatype_", data.type, ".txt", sep=""),
			row.names = TRUE,
			col.names = NA,
			sep = "\t"
			);

		con <- gzfile(
			paste(out.dir, "/coxph_nodes__", all.training.names, "__datatype_", data.type, "_models.gz", sep=""), 
			"wb"
			);
		serialize(coef.nodes.edges[[data.type]][["cox.uv"]], con);
		close(con);
		}

	# lets compute subnetwork feature scores
	all.adjacency.matrices <- get.adjacency.matrix(subnets.file);

	# add "_at" to entrez ids
	# PCB: might want to consider making the whole package independent of this
	for (subnet in names(all.adjacency.matrices)) {
		rownames(all.adjacency.matrices[[subnet]]) <- paste(rownames(all.adjacency.matrices[[subnet]]), "_at", sep="");
		colnames(all.adjacency.matrices[[subnet]]) <- paste(colnames(all.adjacency.matrices[[subnet]]), "_at", sep="");
		}

	# lets compute feature scores
	for (p.threshold in feature.selection.p.thresholds) {
		for (subnet in names(all.adjacency.matrices)) {

			# initialise the scores for each of the model
			# PCB: consider giving the models names (nodes-only, edges-only, nodes+edges) instead of "1", "2", and "3" 
			subnet.scores[["1"]][[subnet]] <- 0;
			subnet.scores[["2"]][[subnet]] <- 0;
			subnet.scores[["3"]][[subnet]] <- 0;

			# go over all data types and sum the scores
			for (data.type in data.types) {
				coxph.nodes <- coef.nodes.edges[[data.type]][["nodes.coef"]];
				coxph.edges.coef <- coef.nodes.edges[[data.type]][["edges.coef"]];
				coxph.edges.P <- coef.nodes.edges[[data.type]][["edges.P"]];
				adjacency.matrix <- all.adjacency.matrices[[subnet]];

				nodes <- rownames(adjacency.matrix);
				nodes <- intersect(
					nodes,
					rownames(coxph.nodes)
					);
				adjacency.matrix <- adjacency.matrix[nodes, nodes];

				# nodes score
				# NOTE: for ordinal data types, feature selection is conducted on the 
				# beta/coef for smallest P value group/level
				nodes.valid <- which( coxph.nodes[nodes, "P"] <= p.threshold );
				if (length(nodes.valid) > 0) {
					nodes.valid <- nodes[nodes.valid];
					nodes.coef <- coxph.nodes[nodes.valid, "coef"];
					nodes.coef <- sum(abs(nodes.coef));
					subnet.scores[["1"]][[subnet]] <- subnet.scores[["1"]][[subnet]] + nodes.coef;
					subnet.scores[["2"]][[subnet]] <- subnet.scores[["2"]][[subnet]] + nodes.coef;
					}

				# edges score
				if (length(nodes) > 1) { # ignore self interaction matrix
					for (i in 2:nrow(adjacency.matrix)) {
						j.max <- i - 1;
						edge.cols <- which(adjacency.matrix[i, 1:j.max] == 1);
						if (length(edge.cols) > 0) {
							g1 <- nodes[i];
							for (nodes.i in 1:length(edge.cols)) {
								g2 <- nodes[edge.cols[nodes.i]];
								if (!is.na(coxph.edges.P[g1, g2]) & coxph.edges.P[g1, g2] <= p.threshold) {
									edges.coef <- abs(coxph.edges.coef[g1, g2]);
									subnet.scores[["1"]][[subnet]] <- subnet.scores[["1"]][[subnet]] + edges.coef;
									subnet.scores[["3"]][[subnet]] <- subnet.scores[["3"]][[subnet]] + edges.coef;
									}
								}
							}
						}
					}
				}
			}

		# lets write subnetwork feature scores to file system
		for (model in names(subnet.scores)) {
			x <- as.matrix(subnet.scores[[model]]);
			x <- as.matrix( x[order(x[, 1], decreasing = TRUE),] );
			colnames(x) <- c("score");
			write.table(
				x = x,
				file = paste(out.dir, "/top_subnets_score__TRAINING_", all.training.names,"__model_", model, "__PV_", p.threshold, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep="\t"
				);
			}
		}
	}
