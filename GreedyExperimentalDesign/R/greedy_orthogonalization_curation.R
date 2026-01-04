#' Curate More Orthogonal Vectors Greedily
#' 
#' This function takes a set of allocation vectors and pares them down one-by-one
#' by eliminating the vector that can result in the largest reduction in Avg[ |r_ij| ].
#' It is recommended to begin with a set of unmirrored vectors for speed. Then add the mirrors later
#' for whichever subset you wish.
#' 
#' @param W 		A matrix in in the set \eqn{{-1, 1}^{R x n}} which have R allocation vectors for an experiment of sample size n.
#' @param Rmin 		The minimum number of vectors to consider in a design. The default is the true bottom, two.
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
#' res = greedy_orthogonalization_curation(W, Rmin = 3, verbose = FALSE)
#' res$avg_abs_rij_by_R
#' }
#' @export
greedy_orthogonalization_curation = function(W, Rmin = 2, verbose = FALSE){
	assertClass(W, "matrix")
	assertCount(Rmin)
	R = nrow(W)
	assertNumeric(Rmin, upper = R)
	
	if (R > 1000){
		warning("If number of designs > 1000, it goes pretty slow.\n")
	}
	
	n = ncol(W)
	Rij = abs(W %*% t(W)) / n
	diag(Rij) = 0
	res = greedy_orthogonalization_eliminate_cpp(Rij, Rmin, TRUE)
	avg_abs_rijss = res$avg_abs_rijss
	eliminations = res$eliminated
	remaining = res$remaining

	if (verbose){
		cat("iter", 0, "avg abs rijs", round(avg_abs_rijss[1], 4), "\n")
		if (length(eliminations) > 0){
			for (iter_num in 1 : length(eliminations)){
				cat("iter", iter_num, "obs eliminated:", eliminations[iter_num],
					"avg abs rijs", round(avg_abs_rijss[iter_num + 1], 4), "\n")
			}
		}
	}
	#now return the original allocations in order of elimination as well as the rho reductions
	list(
		avg_abs_rij_by_R =
			data.frame(R = R - (0 : (length(avg_abs_rijss) - 1)), avg_abs_rijs = avg_abs_rijss),
		Wsorted =
			W[c(remaining, rev(eliminations)), ]
	)
}
