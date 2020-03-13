#' Begin Gurobi Optimized Search
#' 
#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through $1_{T}$ space for forced balance designs. Make sure you setup Gurobi properly first. This means
#' applying for a license, downloading, installing, registering it on your computer using the \code{grbgetkey} command
#' with the license file in the default directory. Then, in R, add it to the path within R using something like 
#' \code{.jaddClassPath("/gurobi801/win64/lib/gurobi.jar")}. If this is not setup properly - you are in for a world of 
#' pain with cryptic errors! Currently, this method does not return multiple vectors. This will be improved in a later 
#' version. If you want this functionality now, use the hacked-up method \code{gurobi_multiple_designs}.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param objective			The objective function to use when searching design space. This is a string
#' 							with valid values "\code{mahal_dist}" (the default) or "\code{kernel}".
#' @param Kgram				If the \code{objective = kernel}, this argument is required to be an \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @param time_limit_min	The maximum amount of time the optimizer can run for in minutes. The default is \code{5}.
#' @param max_solutions		The maximum number of solutions Gurobi should retain incidentally while looking for one optimal 
#' 							(if possible given the time limit and constraint of the node limit). The default is 
#' 							\code{NULL} for only the best.
#' @param node_limit		The maximum number of nodes Gurobi should explore. Default is \code{NULL} for no limit.
#' @param verbose			Should Gurobi print its log to screen? Default is \code{TRUE}.
#' @param log_file			Log filename for Gurobi e.g. \code{my_log.txt}. Default is \code{""} for no file log. 
#' @return					An object of type \code{optimal_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner and Bracha Blau
#' @export
initGurobiNumericalOptimizationExperimentalDesignObject = function(
		X = NULL, 
		objective = "mahal_dist",
		Kgram = NULL,
		num_cores = 1, 
		time_limit_min = 5, 
		node_limit = NULL, 
		max_solutions = NULL,
		verbose = TRUE,
		log_file = ""){
	
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
	
	verify_objective_function(objective, Kgram, n)
	
	if (!is.null(Kgram)){
		n = nrow(Kgram)
		p = NA
	} else {
		n = nrow(X)
		p = ncol(X)
	}
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	
	if (objective == "abs_sum_diff"){
		stop("The \"abs_sum_diff\" objective type does not work with Gurobi.")
	}
	SinvX = NULL
	
	if (objective == "mahal_dist"){
		if (p < n){
			SinvX = solve(var(X))
		}
	}
	
	#we are about to construct a GurobiNumericalOptimizeExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	error_obj = NULL
	java_obj = .jnew("GurobiNumericalOptimizeExperimentalDesign.GurobiNumericalOptimizeExperimentalDesign")
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setN", as.integer(n))
	if (objective != "kernel"){
		p = ncol(X)
		.jcall(java_obj, "V", "setP", as.integer(p))
	}
#	cat("time limit min: ", as.numeric(time_limit_min), "\n")
	.jcall(java_obj, "V", "setObjective", objective)
	.jcall(java_obj, "V", "setTimeLimitMin", as.numeric(time_limit_min))
	if (!is.null(node_limit)){
		if (node_limit <= 1){
			stop("Node limit must be one or more.")
		}
		.jcall(java_obj, "V", "setNodeLimit", as.numeric(round(node_limit)))
	}	
	if (!is.null(max_solutions)){
		if (max_solutions < 1){
			stop("Max solutions must be one or more.")
		}
		.jcall(java_obj, "V", "setMaxSolutions", as.integer(round(max_solutions)))
	}	
	
	if (!verbose){
		.jcall(java_obj, "V", "turnGurobiLogOff")
	}
	.jcall(java_obj, "V", "setLogFilename", log_file)
	
	#feed in the gram matrix if applicable
	if (!is.null(Kgram)){
		setGramMatrix(java_obj, Kgram)
	} else {
		#feed in the raw data
		for (i in 1 : n){	
			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
		}
		
		#feed in the inverse var-cov matrix
		for (j in 1 : p){
			.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
		}
	}
	
	#now return information as an object (just a list)
	gurobi_numerical_optimization_experimental_design_search = list()
	gurobi_numerical_optimization_experimental_design_search$X = X
	gurobi_numerical_optimization_experimental_design_search$SinvX = SinvX
	gurobi_numerical_optimization_experimental_design_search$n = n
	gurobi_numerical_optimization_experimental_design_search$p = p
	gurobi_numerical_optimization_experimental_design_search$Kgram = Kgram
	gurobi_numerical_optimization_experimental_design_search$objective = objective
	gurobi_numerical_optimization_experimental_design_search$max_solutions = max_solutions
	gurobi_numerical_optimization_experimental_design_search$time_limit_min = time_limit_min
	gurobi_numerical_optimization_experimental_design_search$java_obj = java_obj
	class(gurobi_numerical_optimization_experimental_design_search) = "gurobi_numerical_optimization_experimental_design_search"
	
	#now start search
	startSearch(gurobi_numerical_optimization_experimental_design_search)
	
	#return the final object
	gurobi_numerical_optimization_experimental_design_search
}




#' Find multiple designs using Gurobi
#' 
#' This method searches through $1_{T}$ space using Gurobi's optimization many times.
#' It finds many different solutions by permuting the rows of the design matrix and 
#' rerunning the optimization.
#' 
#' @param X 		The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 					(one for each measurement on the subject). This is the design matrix you wish 
#' 					to search for a more optimal design.
#' @param r 		The number of vectors that should be returned
#' @param ... 		Additional arguments to be passed to \code{initGurobiNumericalOptimizationExperimentalDesignObject}.
#' @return			A matrix of allocation vectors of dimension \code{r x n}.	 
#' 
#' @author Kapelner
#' @export
gurobi_multiple_designs = function(X, r, ...){
	n = nrow(X)
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	for (nsim in 1 : r){		
		random_indices = sample(1 : n)
		X_randomized = X[random_indices, , drop = FALSE]
		gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(X_randomized, ...)
		indicT = resultsGurobiNumericalOptimizeExperimentalDesign(gnoed)$indicT
		indicTs[nsim, ] = indicT[order(random_indices)]
	}	
	indicTs
}

#' Find multiple designs using Gurobi and returns the minimum
#' 
#' This method searches through $1_{T}$ space using Gurobi's optimization many times.
#' It finds many different solutions by permuting the rows of the design matrix and 
#' rerunning the optimization.
#' 
#' @param X 		The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 					(one for each measurement on the subject). This is the design matrix you wish 
#' 					to search for a more optimal design.
#' @param r 		The number of vectors that should be returned
#' @param objective The objective function to calculate. Default is \code{"mahal_dist"} for Mahalanobis distance
#' @param ... 		Additional arguments to be passed to \code{initGurobiNumericalOptimizationExperimentalDesignObject}.
#' @return			A list with the minimum objective value and its vector	 
#' 
#' @author Kapelner
#' @export
gurobi_min_of_multiple_designs = function(X, r, objective = "mahal_dist", ...){
	indicTs = gurobi_multiple_designs(X, r, ...)
	obj_min = .Machine$double.xmax
	i_min = NA
	for (i in 1 : r){
		obj_i = compute_objective_val(X, indicTs[i, ], objective = "mahal_dist")
		if (obj_i < obj_min){
			obj_min = obj_i
			i_min = i
		}
	}
	list(indicT = indicTs[i_min, ], obj = obj_min)
}





#' Query the Gurobi Results
#' 
#' Returns the results (thus far) of the Gurobi numerical optimization design search
#' 
#' @param obj 				The \code{gurobi_numerical_optimization_experimental_design_search} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
resultsGurobiNumericalOptimizeExperimentalDesign = function(obj){
	indicTs = .jcall(obj$java_obj, "[[I", "getIndicTs", simplify = TRUE)
	indicTs = t(unique(indicTs)) #remove all duplicates
	###hack.... some Gurobi solutions are illegal because they do not respect n_T - n_C. Manually remove these
	indicTs = indicTs[, colSums(indicTs) == obj$n / 2]
	if (is.null(obj$max_solutions)){ #we only wanted one

	} else {
		
	}
	obj_vals = apply(indicTs, 2, FUN = function(wj){compute_objective_val(X, wj, objective = obj$objective)})
	indicTs = indicTs[, order(obj_vals), drop = FALSE]
	list(indicTs = indicTs, obj_vals = sort(obj_vals))	
}




