calculate.network.coefficients <- function(data.directory = ".", output.directory = ".", training.datasets = NULL, data.types = c("mRNA"), subnets.file.flattened = NULL, subset = NULL) {

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

	cancer.data <- load.cancer.datasets(
		datasets.to.load = training.datasets, 
		data.types = data.types, 
		data.directory = data.directory,
		subset = subset
		);

	# DATA PROCESSING
	scaled.data <- cancer.data;

	# SCALE DATA
	for (data.type in data.types) {
		for (i in 1:length(scaled.data[["all.data"]][[data.type]]) ) {
			scaled.data[["all.data"]][[data.type]][[i]] <- as.data.frame( t( scale( t(scaled.data[["all.data"]][[data.type]][[i]]) ) ) );
			}
		}

	# ANALYZE EACH SAMPLE
	# make the ProbeSet IDs
	gene.pairs$ProbeID1 <- paste(as.character(gene.pairs$GeneID1), "_at", sep = "");
	gene.pairs$ProbeID2 <- paste(as.character(gene.pairs$GeneID2), "_at", sep = "");

	# initialise return oject
	nodes.edges.stats <- list();

	for (data.type in data.types) {

		coxph.nodes <- matrix(
			data = NA,
			nrow = length(subnets.genes),
			ncol = 2,
			dimnames = list(
				subnets.genes,
				c("coef", "P")
				)
			);

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

			results <- fit.interaction.model(
				feature1 = gene.pairs$ProbeID1[i],
				feature2 = gene.pairs$ProbeID2[i],
				expression.data = scaled.data[["all.data"]][[data.type]],
				survival.data = scaled.data$all.survobj
				);

			# sometimes same gene is in multiple interactions and hence return all NULL row
			# and with other interactions, return numeric HR and we would like to keep numeric
			# HRs not NA when numeric HR is available
			if (!is.na(results[1])) { 
				coxph.nodes[gene.pairs$ProbeID1[i], "coef"] <- log2(results[1]);
				coxph.nodes[gene.pairs$ProbeID1[i], "P"]  <- results[2];
				coxph.nodes[gene.pairs$ProbeID2[i], "coef"] <- log2(results[3]);
				coxph.nodes[gene.pairs$ProbeID2[i], "P"]  <- results[4];
				}

			# store the edges HR and P in a matrix - for faster lookups
			coxph.edges.coef[gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <- log2(results[5]);
			coxph.edges.coef[gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <- log2(results[5]);
			coxph.edges.P[ gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <- results[6];
			coxph.edges.P[ gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <- results[6];

			}

		# populate return object
		nodes.edges.stats[[data.type]][["nodes.coef"]] <- coxph.nodes;
		nodes.edges.stats[[data.type]][["edges.coef"]] <- coxph.edges.coef;
		nodes.edges.stats[[data.type]][["edges.P"]] <- coxph.edges.P;

		}

	return(nodes.edges.stats);

	}
