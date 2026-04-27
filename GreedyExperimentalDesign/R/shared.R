#' Starts the parallelized greedy design search. 
#' 
#' Once begun, this function cannot be run again.
#' 
#' @param obj 		The \code{experimental_design} object that will be running the search
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 3,
#'   num_cores = 1,
#'   start = FALSE,
#'   wait = FALSE,
#'   objective = "abs_sum_diff",
#'   verbose = FALSE
#' )
#' startSearch(ged)
#' stopSearch(ged)
#' }
#' @export
startSearch = function(obj){
	if (.jcall(obj$java_obj, "Z", "began")){
		stop("Search Already begun.")
	}
	.jcall(obj$java_obj, "V", "beginSearch")
}

#' Stops the parallelized greedy design search. 
#' 
#' Once stopped, it cannot be restarted.
#' 
#' @param obj 		The \code{experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 3,
#'   num_cores = 1,
#'   start = TRUE,
#'   wait = FALSE,
#'   objective = "abs_sum_diff",
#'   verbose = FALSE
#' )
#' stopSearch(ged)
#' }
#' @export
stopSearch = function(obj){
	.jcall(obj$java_obj, "V", "stopSearch")
}

#' Generates a design matrix with standardized predictors. 
#' 
#' This function is useful for debugging.
#' 
#' @param n					Number of rows in the design matrix 
#' @param p 				Number of columns in the design matrix
#' @param covariate_gen		The function to use to draw the covariate realizations (assumed to be iid).
#' 							This defaults to \code{rnorm} for $N(0,1)$ draws.
#' @param ...				Optional arguments to be passed to the \code{covariate_dist} function.
#' @return 					THe design matrix
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' X = generate_stdzied_design_matrix(n = 6, p = 2)
#' colMeans(X)
#' }
#' @export
generate_stdzied_design_matrix = function(n = 50, p = 1, covariate_gen = rnorm, ...){
	X = matrix(covariate_gen(n * p, ...), nrow = n, ncol = p)
	#now standardize the matrix to make things easier later
	standardize_data_matrix(X)
}


#' Returns the amount of time elapsed
#' 
#' @param obj 		The \code{experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 1,
#'   num_cores = 1,
#'   start = TRUE,
#'   wait = TRUE,
#'   objective = "abs_sum_diff",
#'   verbose = FALSE
#' )
#' searchTimeElapsed(ged)
#' }
#' @export
searchTimeElapsed = function(obj){
	.jcall(obj$java_obj, "I", "timeElapsedInSeconds")
}

#' Computes a numerically stable inverse of a covariance matrix
#'
#' @param X			The n x p design matrix
#' @param ridge		Initial ridge penalty added to the diagonal
#' @param max_ridge_steps	Maximum number of ridge escalation attempts
#' @return			The inverse covariance matrix
#'
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' X = matrix(rnorm(20), nrow = 10)
#' Sinv = safe_cov_inverse(X)
#' dim(Sinv)
#' }
#' @export
safe_cov_inverse = function(X, ridge = 1e-8, max_ridge_steps = 6){
	S = stats::var(X)
	if (!all(is.finite(S))){
		stop("Covariance matrix contains non-finite values.")
	}
	Sinv = tryCatch(solve(S), error = function(e) NULL)
	if (!is.null(Sinv)){
		return(Sinv)
	}
	ridge_val = ridge
	for (i in seq_len(max_ridge_steps)){
		Sinv = tryCatch(solve(S + diag(ridge_val, ncol(S))), error = function(e) NULL)
		if (!is.null(Sinv)){
			return(Sinv)
		}
		ridge_val = ridge_val * 10
	}
	sv = svd(S)
	if (length(sv$d) == 0){
		stop("Covariance matrix is empty.")
	}
	tol = max(dim(S)) * max(sv$d) * .Machine$double.eps
	d_inv = ifelse(sv$d > tol, 1 / sv$d, 0)
	sv$v %*% (d_inv * t(sv$u))
}

#' Computes Objective Value From Allocation Vector
#' 
#' Returns the objective value given a design vector as well an an objective function.
#' This is sometimes duplicated in Java. However, within Java, tricks are played to make
#' optimization go faster so Java's objective values may not always be the same as the true
#' objective function (e.g. logs or constants dropped).
#' 
#' @param X 		 	The n x p design matrix
#' @param indic_T		The n-length binary allocation vector
#' @param objective		The objective function to use. Default is \code{abs_sum_diff} and the other option is 
#' 						\code{mahal_dist}.
#' @param inv_cov_X		Optional: the inverse sample variance covariance matrix. Use this
#' 						argument if you will be doing many calculations since passing this
#' 						in will cache this data.
#' @param use_safe_inverse	Should a regularized inverse be used for the Mahalanobis objective?
#' 							Default is \code{FALSE}.
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' X = matrix(rnorm(12), nrow = 6)
#' indic_T = c(1, 0, 1, 0, 1, 0)
#' compute_objective_val(X, indic_T, objective = "abs_sum_diff")
#' }
#' @export
compute_objective_val = function(X, indic_T, objective = "abs_sum_diff", inv_cov_X = NULL, use_safe_inverse = FALSE){
	assertLogical(use_safe_inverse)
	if (!isTRUE(all.equal(sort(unique(indic_T)), c(0, 1)))){
		stop("indic_T must be binary")
	}
	X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_C = X[indic_T == 0, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_T_bar = colMeans(X_T)
	X_C_bar = colMeans(X_C)	
	
	if (objective == "abs_sum_diff"){
		s_j_s = apply(X, 2, sd)
		sum(abs((X_T_bar - X_C_bar) / s_j_s))
	} else if (objective == "mahal_dist"){
		#saves computation to pass it in if you're doing a lot of them in a row
		if (is.null(inv_cov_X)){
			if (use_safe_inverse){
				inv_cov_X = safe_cov_inverse(X)
			} else {
				inv_cov_X = solve(stats::var(X))
			}
		}	
		X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar) #need to matricize for next computation
		as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)
	} else {
		stop("objective invalid.")
	}
}

#' Standardizes the columns of a data matrix.
#' 
#' @param X 		 	The n x p design matrix
#' @return				The n x p design matrix with columns standardized
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' X = matrix(rnorm(12), nrow = 6)
#' Xstd = standardize_data_matrix(X)
#' colMeans(Xstd)
#' }
#' @export
standardize_data_matrix = function(X){
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
}

#private
verify_objective_function = function(objective, Kgram = NULL, n = NULL){
	if (objective != "mahal_dist" && objective != "abs_sum_diff" && objective != "kernel"){
		stop("Objective function must be one of the following:\n  mahal_dist\n  abs_sum_diff\n  kernel\n\n")
	}
	if (objective == "kernel"){
		if (is.null(Kgram) || is.null(n)){
			stop("You must specify a gram matrix \"Kgram\" and \"n\".\n")
		}
		if (!inherits(Kgram, "kernelMatrix") && !inherits(Kgram, "matrix")){
			stop("The gram matrix must be type kernelMatrix or type matrix.\n")
		}
		if (!all.equal(dim(Kgram), c(n, n))){
			stop("The gram matrix must have dimension n x n.\n")
		}
	}
	if (!is.null(Kgram) && objective != "kernel"){
		stop("If you specify a gram matrix, you must specify the \"kernel\" objective.\n")
	}
}




####################C++ functions

#' Generates homogeneous block design allocations rapidly
#' 
#' @param B 		 	The number of blocks in the design
#' @param nR			The number of allocation vectors
#' @param dummy_block   The subvector of allocations in each block that will be permuted
#' @return				A matrix with rows being the nR random block allocation of sample size B x length(dummy_block).
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' generate_block_design_cpp_wrap(B = 2, nR = 3, dummy_block = c(1, 0))
#' }
#' @export
generate_block_design_cpp_wrap <- function(B, nR, dummy_block) {
    do.call(cbind, .Call('_GreedyExperimentalDesign_generate_block_design_cpp', PACKAGE = 'GreedyExperimentalDesign', B, nR, dummy_block))
}


#' Computes a Euclidean-squared distance matrix rapidly
#' 
#' @param X 		 	A numeric matrix with n rows representing each subject and p columns which are measurements on each subject
#' @return				The n x n Euclidean distances squared
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' X = matrix(c(0, 1, 2, 3), nrow = 2)
#' compute_distance_matrix_cpp_wrap(X)
#' }
#' @export
compute_distance_matrix_cpp_wrap <- function(X) {
    .Call('_GreedyExperimentalDesign_compute_distance_matrix_cpp', PACKAGE = 'GreedyExperimentalDesign', X)
}

#' Shuffles a vector rapidly
#' 
#' @param w 		 	The vector to be shuffled
#' @param seed		Optional integer seed; use NA to draw from the system clock
#' @return				The vector with elements shuffled
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' shuffle_cpp_wrap(1:5, seed = 1)
#' }
#' @export
shuffle_cpp_wrap <- function(w, seed = NA_integer_) {
    .Call('_GreedyExperimentalDesign_shuffle_cpp', PACKAGE = 'GreedyExperimentalDesign', w, seed)
}

#' Tests if a vector has all elements the same
#' 
#' @param w 		 	The vector to be queried
#' @return				A boolean if it has all same elements
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' all_elements_same_cpp_wrap(c(1, 1, 1))
#' all_elements_same_cpp_wrap(c(1, 2, 1))
#' }
#' @export
all_elements_same_cpp_wrap <- function(w) {
    .Call('_GreedyExperimentalDesign_all_elements_same_cpp', PACKAGE = 'GreedyExperimentalDesign', w)
}

#' Create all binary Y's convenience function using a randomized design
#' 
#' @param pCs			Control-group success probabilities (length \code{two_n})
#' @param pTs			Treatment-group success probabilities (length \code{two_n})
#' @param W				Assignment matrix with \code{nY} rows and \code{two_n} columns
#' @param two_n			Total number of units
#' @param nY				Number of Y vectors to generate
#' @return				A matrix of boolean Y's
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' pCs = rep(0.2, 4)
#' pTs = rep(0.8, 4)
#' W = matrix(c(1, 0, 1, 0, 0, 1, 0, 1), nrow = 2, byrow = TRUE)
#' create_all_ys_cpp_wrap(pCs, pTs, W, two_n = 4, nY = 2)
#' }
#' @export
create_all_ys_cpp_wrap <- function(pCs, pTs, W, two_n, nY) {
    .Call('_GreedyExperimentalDesign_create_all_ys_cpp', PACKAGE = 'GreedyExperimentalDesign', pCs, pTs, W, two_n, nY)
}

#' Create PM designs
#' 
#' @param indicies_pairs	A matrix of n x 2 indicies where each row is a pair of subjects' indicies
#' @param n					Half the number of subjects i.e. the number of pairs
#' @param r					The number of assignments to generate
#' @return					A matrix of r x 2n PM designs of +1/-1 assignments
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' indicies_pairs = matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
#' gen_pm_designs_cpp_wrap(indicies_pairs, n = 2, r = 3)
#' }
#' @export
gen_pm_designs_cpp_wrap <- function(indicies_pairs, n, r) {
    .Call('_GreedyExperimentalDesign_gen_pm_designs_cpp', PACKAGE = 'GreedyExperimentalDesign', indicies_pairs, n, r)
}

# internal helper to avoid hard failure when older Java jars omit setVerbose
set_verbose_if_available <- function(java_obj, verbose) {
	tryCatch(
		.jcall(java_obj, "V", "setVerbose", verbose),
		error = function(e) {
			if (grepl("setVerbose", conditionMessage(e), fixed = TRUE)) {
				return(invisible(FALSE))
			}
			stop(e)
		}
	)
	invisible(TRUE)
}

# internal GPU dispatch flag
ged_gpu_option <- function() {
	isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))
}

ged_use_gpu <- function(require_available = FALSE) {
	use_gpu = ged_gpu_option()
	if (use_gpu && require_available && !isTRUE(ged_gpu_available())) {
		warning("GreedyExperimentalDesign.use_gpu is TRUE, but no GPU backend is available. Falling back to CPU.")

		return(FALSE)
	}
	use_gpu
}

# internal helper to avoid hard failure when older Java jars omit setUseGpu.
# The raw option flag is passed to Java; R-side dispatch separately checks
# native backend availability before using GPU paths.
set_use_gpu_if_available <- function(java_obj, use_gpu = ged_gpu_option()) {
	tryCatch(
		.jcall(java_obj, "V", "setUseGpu", use_gpu),
		error = function(e) {
			if (grepl("setUseGpu", conditionMessage(e), fixed = TRUE)) {
				return(invisible(FALSE))
			}
			stop(e)
		}
	)
	invisible(TRUE)
}

# internal helper to avoid unnecessary coercion for w_diff
compute_indicTs_from_pairs <- function(pairs, w_diff, n) {
	if (is.logical(w_diff)) {
		compute_indicTs_from_pairs_lgl_cpp(pairs, w_diff, n)
	} else if (is.raw(w_diff)) {
		compute_indicTs_from_pairs_raw_cpp(pairs, w_diff, n)
	} else {
		if (!anyNA(w_diff) && all(w_diff == 0L | w_diff == 1L)) {
			return(compute_indicTs_from_pairs_lgl_cpp(pairs, w_diff == 1L, n))
		}
		compute_indicTs_from_pairs_cpp(pairs, w_diff, n)
	}
}
