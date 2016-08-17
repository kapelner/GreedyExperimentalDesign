#' This method creates an object of type lPSolve_experimental_design and will immediately initiate
#' a search through $1_{T}$ space.
#' 
#' @param X							The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 									(one for each measurement on the subject). This is the design matrix you wish 
#' 									to search for a more optimal design.
#' @param obj_val_cutoff_to_include	Only allocation vectors with objective values lower than this threshold will be returned.
#' 									The default is \code{NULL} which means all vectors are returned.
#' @param max_designs 				The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 									so you can search however long you wish as the search can be stopped at any time by
#' 									using the \code{\link{stopSearch}} method 
#' @param objective					The objective function to use when searching the design space. This is a string
#' 									"\code{abs_sum_diff}" (default) or "\code{mahal_dist}."
#' @param wait						Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 									default is \code{FALSE}.
#' @param start						Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores 				The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return							An object of type \code{lPSolve_experimental_design_search} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
initLPSolveExperimentalDesignObject = function(X, 
		objective = "mahal_dist", 
		wait = FALSE, 
		start = TRUE,
		num_cores = 1){
	
	#get dimensions immediately
	n = nrow(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	p = ncol(X)
	
	if (objective == "abs_sum_diff"){
		#standardize it -- much faster here
		Xstd = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
	}
	if (objective == "mahal_dist"){
		if (p < n){
			SinvX = solve(var(X))
		}
	}
	
	#we are about to construct a LPSolveExperimentalDesign java object. First, let R garbage collect
	#to clean up previous LPSolveExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("LPSolveExperimentalDesign.LPSolveExperimentalDesign")
	.jcall(java_obj, "V", "setMaxDesigns", as.integer(max_designs))
	if (!is.null(obj_val_cutoff_to_include)){
		.jcall(java_obj, "V", "setObjValCutoffToInclude", obj_val_cutoff_to_include)
	}
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setNandP", as.integer(n), as.integer(p))
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}	
	
	#feed in the data
	for (i in 1 : n){	
		if (objective == "abs_sum_diff"){
			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
		} else {
			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
		}
	}
	
	#feed in the inverse var-cov matrix
	if (objective == "mahal_dist"){
		if (p < n){
			for (j in 1 : p){
				.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
			}
		}
	}
	
	#now return information as an object (just a list)
	lPSolve_experimental_design_search = list()
	lPSolve_experimental_design_search$max_designs = max_designs
	lPSolve_experimental_design_search$obj_val_cutoff_to_include = obj_val_cutoff_to_include
	lPSolve_experimental_design_search$start = start
	lPSolve_experimental_design_search$wait = wait
	lPSolve_experimental_design_search$X = X
	lPSolve_experimental_design_search$n = n
	lPSolve_experimental_design_search$p = p
	lPSolve_experimental_design_search$objective = objective
	lPSolve_experimental_design_search$java_obj = java_obj
	class(lPSolve_experimental_design_search) = "lPSolve_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(lPSolve_experimental_design_search)
	}
	#return the final object
	lPSolve_experimental_design_search
}

#' Prints a summary of a \code{lPSolve_experimental_design_search} object
#' 
#' @param x			The \code{lPSolve_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print lPSolve_experimental_design_search
#' @export
print.lPSolve_experimental_design_search = function(x, ...){
	progress = lPSolveSearchCurrentProgress(x)
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the LPSolveExperimentalDesign. Did you run \"startSearch?\"\n")
	} else if (progress == x$max_designs){
		cat("The search completed in", time_elapsed, "seconds.", progress, "vectors have been found.\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%) in ", time_elapsed," seconds.\n", sep = "")
	}
}

#' Prints a summary of a \code{lPSolve_experimental_design_search} object
#' 
#' @param object		The \code{lPSolve_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary lPSolve_experimental_design_search
#' @export
summary.lPSolve_experimental_design_search = function(object, ...){
	print(object, ...)
}

# Returns the number of vectors found by the lPSolve design search
# 
# @param obj 		The \code{lPSolve_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
lPSolveSearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}


#' Returns the results (thus far) of the lPSolve design search
#' 
#' @param obj 			The \code{lPSolve_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
resultsLPSolveSearch = function(obj){
	obj_vals = NULL
	if (!is.null(obj$obj_val_cutoff_to_include)){
		obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	}	
	ending_indicTs = sapply(.jcall(obj$java_obj, "[[I", "getEndingIndicTs"), .jevalArray)
	
	lPSolve_experimental_design_search_results = list(
			obj_vals = obj_vals, 
			ending_indicTs = ending_indicTs
	)
	class(lPSolve_experimental_design_search_results) = "lPSolve_experimental_design_search_results"
	#return the final object
	lPSolve_experimental_design_search_results
}
#}