#' Computes Randomization Metrics (explained in paper) about a design algorithm
#' 
#' @param designs	A matrix where each column is one design.
#' 
#' @return 			A list of resulting data: the probability estimates for
#' 					each pair in the design of randomness where estmates close
#' 					to ~0.5 represent random assignment, then the entropy metric
#' 					and the distance metric.
#' 
#' @author Adam Kapelner
#' @export
compute_randomization_metrics = function(designs){
	n = nrow(designs)
	r = ncol(designs)
	
	gc() #Delete at your own risk!
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("DesignMetrics.RandomizationMetrics")
	.jcall(java_obj, "V", "setNandR", as.integer(n), as.integer(r))
#	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	
	#feed in the data
	for (j in 1 : r){
		.jcall(java_obj, "V", "setDesign", as.integer(j - 1), as.integer(designs[, j])) #java indexes from 0...n-1
	}
	#get it going
	.jcall(java_obj, "V", "compute")
	
	#harvest the data and return it as a list
	p_hat_ijs = sapply(.jcall(java_obj, "[[D", "getPhats"), .jevalArray)
	rand_entropy_metric = .jcall(java_obj, "D", "getRandEntropyMetric")
	rand_norm_se_metric = .jcall(java_obj, "D", "getRandStdErrMetric")
	list(
		p_hat_ijs = p_hat_ijs, 
		rand_entropy_metric = rand_entropy_metric, 
		rand_norm_se_metric = rand_norm_se_metric
	)
}