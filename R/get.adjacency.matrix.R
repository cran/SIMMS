#' A utility function to convert tab delimited networks file into adjacency
#' matrices
#' 
#' A utility function to convert tab-delimited networks file into adjacency
#' matrices
#' 
#' 
#' @param subnets.file A tab-delimited file containing networks. New networks
#' start with a new line with '#' at the begining of network name and
#' subsequent lines contain a binary interaction per line
#' @return A list of adjacency matrices
#' @author Syed Haider
#' @keywords Networks
#' @examples
#' 
#' subnets.file <- get.program.defaults()[["subnets.file"]];
#' all.adjacency.matrices <- get.adjacency.matrix(subnets.file);
#' 
#' @export get.adjacency.matrix
get.adjacency.matrix <- function(subnets.file = NULL) {

	all.adjacency.matrices <- list();
	subnets <- readLines(subnets.file, ok = TRUE);

	graph.name <- "";
	vertices <- "";
	interactions <- "";

	for(i in seq(1,length(subnets),1)) {

		# check if its a header line 
		if (length(grep("^#", subnets[i], perl = TRUE)) > 0) {

			# time to process previous subgraph
			if (nchar(as.character(vertices)) > 0) {

				# make a matrix of this graph
				adjacency.matrix <- make.matrix(vertices, interactions);

				all.adjacency.matrices[[graph.name]] <- adjacency.matrix;

				# reinitialise everything else
				vertices <- "";
				interactions <- "";
				}
			graph.name <- make.names(gsub("\t$", "", subnets[i]));
			}
		else {
			id.p1.p2 <- unlist(strsplit(subnets[i], "\t"));
			id.p1.p2 <- gsub("\\(|\\)", "-", id.p1.p2, perl=TRUE);
			p1 <- paste("\"",id.p1.p2[2],"\"", sep="");
			p2 <- paste("\"",id.p1.p2[3],"\"", sep="");

			# this vertex is not already seen in this sub graph, lets add
			if (length(grep(p1, vertices, perl = TRUE)) < 1) {
				vertices <- paste(vertices, p1, sep=",");
				}
			if (length(grep(p2, vertices, perl = TRUE)) < 1) {
				vertices <- paste(vertices, p2, sep = ",");
				}
			interactions <- paste(interactions, ",\"", id.p1.p2[2],":", id.p1.p2[3],"\"", sep="");
			}

		}

	# make a matrix of this graph
	adjacency.matrix <- make.matrix(vertices, interactions);

	all.adjacency.matrices[[graph.name]] <- adjacency.matrix;
	
	return (all.adjacency.matrices);

	}
