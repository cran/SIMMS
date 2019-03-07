#' A utility function to return the inst/ directory of the installed package
#' and other default settings
#' 
#' A utility function to return the inst/ directory of the installed package to
#' get the test datasets and other program related data contents
#' 
#' 
#' @param networks.database Name of the pathway networks database. Default to
#' NCI PID/Reactome/Biocarta i-e "default"
#' @return Returns a list of paths to the input directories/files where the
#' contents of this package are installed
#' @author Syed Haider
#' @keywords IO
#' @examples
#' 
#' program.data <- get.program.defaults();
#' 
#' @export get.program.defaults
get.program.defaults <- function(networks.database = "default") {

	# template file to search for
	entry.secret = "pathway_based_sub_networks.txt";	

	# make a list of potential locations for the datasets file
	program.data.dirs <- paste(.libPaths(), "/SIMMS/programdata/networkdb/", networks.database, "/", sep = "");

	# then search all locations
	file.checks <- file.exists( paste(program.data.dirs, "/", entry.secret, sep = "") );

	# check to see if the file was actually found
	if (any(file.checks)) {
		program.data.dir <- program.data.dirs[ order(file.checks, decreasing = TRUE)[1] ];
		}
	else {
		stop("Unable to find pathway_based_sub_networks.txt file");
		}
	
	return (
		list(
			"program.data.dir" = program.data.dir,
			"subnets.file" = paste(program.data.dir, "/pathway_based_sub_networks.txt", sep=""),
			"subnets.file.flattened" = paste(program.data.dir, "/pathway_based_networks_flattened.txt", sep=""),
			"subnets.file.all" = paste(program.data.dir, "/pathway_based_sub_networks_all.txt", sep=""), # only used by BL.pipeline.SIMMS
			"test.data.dir" = paste(program.data.dir, "../../testdata/", sep="")
			)
		);
	}
