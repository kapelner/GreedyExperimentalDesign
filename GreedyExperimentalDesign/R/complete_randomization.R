#' Implements complete balanced randomization
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization = function(n, r){
	X = generate_stdzied_design_matrix(n = n, p = 1)
	rd = initRerandomizationExperimentalDesignObject(X, max_designs = r)
	res = resultsRerandomizationSearch(rd)
	res$ending_indicTs
}