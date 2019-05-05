#' Calculate Cox statistics for input dataset
#' 
#' Function to compute hazard ratios for the genes in pathway-derived networks,
#' by aggregating input datasets into one training cohort. The hazard ratios
#' are computed for each pair by calculating the HR of each gene independently
#' and as an interaction (i.e. y = HR(A) + HR(B) + HR(A:B)
#' 
#' 
#' @param data.directory Path to the directory containing datasets as specified
#' by \code{training.datasets}
#' @param output.directory Path to the output folder where intermediate and
#' results files will be saved
#' @param training.datasets A vector containing names of training datasets
#' @param data.types A vector of molecular datatypes to load. Defaults to
#' c('mRNA')
#' @param data.types.ordinal A vector of molecular datatypes to be treated as
#' ordinal. Defaults to c('cna')
#' @param centre.data A character string specifying the centre value to be used for 
#' scaling data. Valid values are: 'median', 'mean', or a user defined numeric threshold
#' e.g. '0.3' when modelling methylation beta values. This value is used for both scaling
#' as well as for dichotomising data for estimating univariate betas from Cox model.
#' Defaults to 'median'
#' @param subnets.file.flattened File containing all the binary ineractions
#' derived from pathway-derived networks
#' @param truncate.survival A numeric value specifying survival truncation in
#' years. Defaults to 100 years which effectively means no truncation
#' @param subset A list with a Field and Entry component specifying a subset of
#' patients to be selected whose annotation Field matches Entry
#' @return Returns a list of matrices for each of the data types. Matrices
#' contain nodes HR/P, edges HR and edges P.
#' @author Syed Haider & Paul C. Boutros
#' @keywords survival
#' @examples
#' 
#' options("warn" = -1);
#' program.data <- get.program.defaults(networks.database = "test");
#' data.directory <- program.data[["test.data.dir"]];
#' subnets.file.flattened <- program.data[["subnets.file.flattened"]];
#' coef.nodes.edges <- calculate.network.coefficients(
#'   data.directory = data.directory,
#'   output.directory = ".",
#'   training.datasets = c("Breastdata1"),
#'   data.types = c("mRNA"),
#'   subnets.file.flattened = subnets.file.flattened
#'   );
#' 
#' @export calculate.network.coefficients
calculate.network.coefficients <- function(data.directory = ".", output.directory = ".", training.datasets = NULL, data.types = c("mRNA"), data.types.ordinal = c("cna"), centre.data = "median", subnets.file.flattened = NULL, truncate.survival = 100, subset = NULL) {

	all.training.names <- paste(sort(training.datasets), collapse="_");
	gene.pairs <- read.table(
		file = subnets.file.flattened,
		header = TRUE,
		row.names = NULL,
		sep = "\t",
		as.is = TRUE
		);

	# focus only on unique records
	gene.pairs <- unique(gene.pairs);
	subnets.genes <- unique(c(gene.pairs$GeneID1, gene.pairs$GeneID2));
	subnets.genes <- paste(subnets.genes, "_at", sep = "");
	# PCB: does the above line this only works with Affymetrix data or data with an _at trailing and a unique ID before it?

	cancer.data <- SIMMS::load.cancer.datasets(
		truncate.survival = truncate.survival,
		datasets.to.load = training.datasets, 
		data.types = data.types, 
		data.directory = data.directory,
		subset = subset
		);

	# ANALYZE EACH SAMPLE
	# make the ProbeSet IDs
	gene.pairs$ProbeID1 <- paste(as.character(gene.pairs$GeneID1), "_at", sep = "");
	gene.pairs$ProbeID2 <- paste(as.character(gene.pairs$GeneID2), "_at", sep = "");

	# initialise return oject
	nodes.edges.stats <- list();

	for (data.type in data.types) {

		data.type.ordinal <- FALSE;
		if (data.type %in% data.types.ordinal) {
			data.type.ordinal <- TRUE;
			}

		coxph.nodes <- matrix(
			data = NA,
			nrow = length(subnets.genes),
			ncol = 2,
			dimnames = list(
				subnets.genes,
				c("coef", "P")
				)
			);

		coxph.nodes.obj <- list();

		coxph.edges.coef <- matrix(
			data = 1,
			nrow = length(subnets.genes),
			ncol = length(subnets.genes),
			dimnames = list(
				subnets.genes,
				subnets.genes
				)
			);

		# just copy the structure to get the p-value matrix
		coxph.edges.P <- coxph.edges.coef;

		# fill in the empty objects
		for (i in 1:nrow(gene.pairs)) {

			results <- SIMMS::fit.interaction.model(
				feature1 = gene.pairs$ProbeID1[i],
				feature2 = gene.pairs$ProbeID2[i],
				expression.data = cancer.data[["all.data"]][[data.type]],
				survival.data = cancer.data$all.survobj,
				data.type.ordinal = data.type.ordinal,
				centre.data = centre.data
				);

			#cat("\nNew pair\t", gene.pairs$ProbeID1[i], "\t", gene.pairs$ProbeID2[i]);
			#print(results);

			# isolate results of g1 and g2
			if (is.list(results)) {

				if (!data.type.ordinal) {
					results[["cox.uv.1"]][["cox.stats"]][1] <- log2(
						results[["cox.uv.1"]][["cox.stats"]][1]
						);
					results[["cox.uv.2"]][["cox.stats"]][1] <- log2(
						results[["cox.uv.2"]][["cox.stats"]][1]
						);
					results.g1 <- results[["cox.uv.1"]][["cox.stats"]];
					results.g2 <- results[["cox.uv.2"]][["cox.stats"]];
					}
				else {
					results.gx <- list();
					for (cox.uv.i in c("cox.uv.1", "cox.uv.2")) {

						# check if cox model failed
						if (length(results[[cox.uv.i]][["cox.obj"]]) > 1 && !is.na(results[[cox.uv.i]][["cox.obj"]])) {
							results[[cox.uv.i]][["cox.stats"]][, "HR"] <- log2(
								results[[cox.uv.i]][["cox.stats"]][, "HR"]
								);
							colnames(results[[cox.uv.i]][["cox.stats"]])[1] <- "coef";

							# if ordinal data type - pick the smallest P
							min.p <- which(
								results[[cox.uv.i]][["cox.stats"]][, "p"] ==
									min(results[[cox.uv.i]][["cox.stats"]][, "p"], na.rm = T)
								);
							if (length(min.p) > 0) {
								results.gx[[cox.uv.i]] <- results[[cox.uv.i]][["cox.stats"]][min.p, ];
								}
							else {
								results.gx[[cox.uv.i]] <- NA;
								}
							}
						else {
							results.gx[[cox.uv.i]] <- NA;
							}
						}
					results.g1 <- results.gx[["cox.uv.1"]];
					results.g2 <- results.gx[["cox.uv.2"]];
					}

				# sometimes same gene is in multiple interactions and hence return all NULL row
				# and with other interactions, return numeric HR and we would like to keep numeric
				# HRs not NA when numeric HR is available
				if (!is.na(results.g1[1])) { 

					# store HR and P
					coxph.nodes[gene.pairs$ProbeID1[i], "coef"] <- results.g1[1];
					coxph.nodes[gene.pairs$ProbeID1[i], "P"]  <- results.g1[4];

					# store cox fit objects
					coxph.nodes.obj[[gene.pairs$ProbeID1[i]]] <- results[["cox.uv.1"]][["cox.stats"]];
					}

				if (!is.na(results.g2[1])) { 

					# store HR and P
					coxph.nodes[gene.pairs$ProbeID2[i], "coef"] <- results.g2[1];
					coxph.nodes[gene.pairs$ProbeID2[i], "P"]  <- results.g2[4];

					# store cox fit objects
					coxph.nodes.obj[[gene.pairs$ProbeID2[i]]] <- results[["cox.uv.2"]][["cox.stats"]];				
					}

				# store the edges HR and P in a matrix - for faster lookups
				coxph.edges.coef[gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <-
					log2(results[["cox.int"]][1]);
				coxph.edges.coef[gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <-
					log2(results[["cox.int"]][1]);
				coxph.edges.P[ gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <- 
					results[["cox.int"]][2];
				coxph.edges.P[ gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <- 
					results[["cox.int"]][2];

				}
			}

		# populate return object
		nodes.edges.stats[[data.type]][["nodes.coef"]] <- coxph.nodes;
		nodes.edges.stats[[data.type]][["edges.coef"]] <- coxph.edges.coef;
		nodes.edges.stats[[data.type]][["edges.P"]] <- coxph.edges.P;
		nodes.edges.stats[[data.type]][["cox.uv"]] <- coxph.nodes.obj;
		}

	return(nodes.edges.stats);

	}
