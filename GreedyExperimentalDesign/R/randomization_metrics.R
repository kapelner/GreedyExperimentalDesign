#' Computes Randomization Metrics (explained in paper) about a design algorithm
#' 
#' @param designs	A matrix where each column is one design.
#' 
#' @return 			A list of resulting data: the probability estimates for
#' 					each pair in the design of randomness where estmates close
#' 					to ~0.5 represent random assignment, then the entropy metric
#' 					the distance metric, the maximum eigenvalue of the allocation
#' 					var-cov matrix (operator norm) and the squared Frobenius norm 
#' 					(the sum of the squared eigenvalues)
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' designs = matrix(c(1, 0, 1, 0, 0, 1, 0, 1), nrow = 4, ncol = 2)
#' compute_randomization_metrics(designs)
#' }
#' @export
compute_randomization_metrics = function(designs){
	n = nrow(designs)
	r = ncol(designs)
	res = compute_randomization_metrics_cpp(designs)
	p_hat_ijs = res$p_hat_ijs
	rand_entropy_metric = res$rand_entropy_metric
	rand_norm_se_metric = res$rand_norm_se_metric
	
	#for the maximum eigenvalue we need to transform the allocation vector to be in {-1, 1}
	designs[designs == 0] = -1
	#now take the eigendecomposition of the variance-covariance matrix of the allocations
	cov_mat = var(t(designs))
	max_eigenval = NA_real_
	frob_norm_sqd = NA_real_
	if (all(is.finite(cov_mat))){
		e_d = eigen(cov_mat, symmetric = TRUE, only.values = TRUE)
		max_eigenval = max(e_d$values)
		frob_norm_sqd = sum(e_d$values^2)
	}
	
	list(
		p_hat_ijs = p_hat_ijs, 
		rand_entropy_metric = rand_entropy_metric, 
		rand_norm_se_metric = rand_norm_se_metric,
		max_eigenval = max_eigenval,
		frob_norm_sqd = frob_norm_sqd
	)
}
