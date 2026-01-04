#' Curate More Orthogonal Vectors Greedily
#' 
#' This function takes a set of allocation vectors and pares them down one-by-one
#' by eliminating the vector that can result in the largest reduction in Avg[ |r_ij| ].
#' It is recommended to begin with a set of unmirrored vectors for speed. Then add the mirrors later
#' for whichever subset you wish.
#' 
#' @param W 		A matrix in \eqn{{-1, 1}^{R x n}} which have R allocation vectors for an experiment of sample size n.
#' @param R0 		The minimum number of vectors to consider in a design. The default is the true bottom, two.
#' @param verbose 	Default is \code{FALSE} but if not, it will print out a message for each iteration.
#' @return 			A list with two elements: (1) \code{avg_abs_rij_by_R} which is a data frame with R - Rmin + 1 rows and 
#' 					columns R and average absolute r_ij and (2) \code{Wsorted} which provides the collection of vectors in
#' 					sorted by best average absolute r_ij in row order from best to worst.
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' W = matrix(sample(c(-1, 1), 6 * 8, replace = TRUE), nrow = 6)
#' W2 = greedy_orthogonalization_curation2(W, R0 = 4, verbose = FALSE)
#' dim(W2)
#' }
#' @export
greedy_orthogonalization_curation2 = function(W, R0 = 100, verbose = FALSE){
	assertClass(W, "matrix")
	assertCount(R0)
	R = nrow(W)
	assertNumeric(R0, upper = R)
	assertLogical(verbose)
	
	Rij = abs(W %*% t(W)) / ncol(W)
	diag(Rij) = 0 #for speed
	res = greedy_orthogonalization_eliminate_cpp(Rij, R0, FALSE)
	if (verbose && length(res$eliminated) > 0){
		for (iter_num in seq_along(res$eliminated)){
			if (iter_num %% 100 == 0){
				cat("iter", iter_num, "\n")
			}
		}
	}
	#now return the subset of original allocations in order of elimination
	W[res$remaining, ]
}
