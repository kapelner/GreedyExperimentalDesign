#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through $1_{T}$ space.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @param time_limit_min	The maximum amount of time the optimizer can run for in minutes. The default is \code{5}.
#' @param node_limt			The maximum number of nodes Gurobi should explore. Default is \code{NULL} for no limit.
#' @return					An object of type \code{optimal_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner and Bracha Blau
#' @export
initGurobiNumericalOptimizationExperimentalDesignObject = function(X, num_cores = 1, time_limit_min = 5, node_limt = NULL){
	
	#we need to check if the user has Gurobi
	gurobi_exists = FALSE
	for (path in .jclassPath()){
		if (length(grep("gurobi.jar", path)) > 0){
			gurobi_exists = TRUE
		}
	}
	
	if (!gurobi_exists){
		stop("You can only use this feature if you have a license for Gurobi\nand the optimizer installed. If so, you should link gurobi.jar via:\n\n .jaddClassPath(\"/<my path>/lib/gurobi.jar\")\n\n and then check to ensure it is properly listed in:\n\n.jclassPath().\n")
	}
	
	#get dimensions immediately
	n = nrow(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	p = ncol(X)
	
	SinvX = solve(var(X))
	
	#we are about to construct a GurobiNumericalOptimizeExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	error_obj = NULL
	java_obj = .jnew("GurobiNumericalOptimizeExperimentalDesign.GurobiNumericalOptimizeExperimentalDesign")
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setNandP", as.integer(n), as.integer(p))
#	cat("time limit min: ", as.numeric(time_limit_min), "\n")
	.jcall(java_obj, "V", "setTimeLimitMin", as.numeric(time_limit_min))
	if (!is.null(node_limt)){
		if (node_limit <= 1){
			stop("Node limit must exceed one.")
		}
		.jcall(java_obj, "V", "setNodeLimit", as.numeric(round(node_limit)))
	}
	
	#feed in the data
	for (i in 1 : n){	
		.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
	}
	
	#feed in the inverse var-cov matrix
	for (j in 1 : p){
		.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
	}
	
		
	#now return information as an object (just a list)
	gurobi_numerical_optimization_experimental_design_search = list()
	gurobi_numerical_optimization_experimental_design_search$X = X
	gurobi_numerical_optimization_experimental_design_search$SinvX = SinvX
	gurobi_numerical_optimization_experimental_design_search$n = n
	gurobi_numerical_optimization_experimental_design_search$p = p
	gurobi_numerical_optimization_experimental_design_search$time_limit_min = time_limit_min
	gurobi_numerical_optimization_experimental_design_search$java_obj = java_obj
	class(gurobi_numerical_optimization_experimental_design_search) = "gurobi_numerical_optimization_experimental_design_search"
	
	#now start search
	startSearch(gurobi_numerical_optimization_experimental_design_search)
	
	#return the final object
	gurobi_numerical_optimization_experimental_design_search
}


#' Returns the results (thus far) of the Gurobi numerical optimization design search
#' 
#' @param obj 				The \code{gurobi_numerical_optimization_experimental_design_search} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
resultsGurobiNumericalOptimizeExperimentalDesign = function(obj){
	indicT = .jcall(obj$java_obj, "[I", "getBestIndicT", .jevalArray)
	list(indicT = indicT)
}
