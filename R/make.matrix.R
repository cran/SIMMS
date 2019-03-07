#' Utility function used by \code{get.adjacency.matrix()}
#' 
#' Utility function used by \code{get.adjacency.matrix()}
#' 
#' 
#' @param vertices Comma separated list of nodes
#' @param interactions Comma separated list of edges
#' @return Returns adjacency matrix
#' @author Syed Haider
#' @keywords Networks
#' @examples
#' 
#' x1 <- make.matrix("a,b,c", "a:b,b:c");
#' 
#' @export make.matrix
make.matrix <- function(vertices, interactions) {

	# remove extra comma at the beginning
	vertices <- sub("^,","",vertices);
	interactions <- sub("^,","",interactions);

	vertices.list <- unlist(strsplit(vertices,","));
	interactions.list <- unlist(strsplit(interactions,","));

	# let retrieve the gene ids i-e vertices names
    captions <- rep("NA", length(vertices.list));
	for (i in seq(1, length(vertices.list), 1)) {
		captions[i] <- gsub("\"", "", vertices.list[i]);
		}

	adjacency.matrix <- matrix(
		data = 0,
		nrow = length(vertices.list),
		ncol = length(vertices.list),
		dimnames = list(
			captions,
			captions
			)
		);

	# lets populate the adjacency matrix
	for (i in 1:length(interactions.list)) {
		interactors <- unlist(strsplit(gsub("\"","", interactions.list[i]), ":"));
		adjacency.matrix[interactors[1], interactors[2]] <- 1;
		adjacency.matrix[interactors[2], interactors[1]] <- 1;
		}

	return(adjacency.matrix);

	}
