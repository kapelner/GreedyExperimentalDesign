#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through $1_{T}$ space.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopOptimalSearch}} method 
#' @param objective			The objective function to use when greedily searching design space. This is a string
#' 							"\code{abs_sum_diff}" (default) or "\code{mahal_dist}."
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return					An object of type \code{optimal_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initOptimalExperimentalDesignObject = function(X,
		objective = "abs_sum_diff", 
		wait = FALSE, 
		start = TRUE,
		num_cores = 1){
	#get dimensions immediately
	n = nrow(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	p = ncol(X)
	
	#standardize it
	Xstd = apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})
	
	if (p <= n){
		SinvXstd = solve(var(Xstd))
	}
	
	#we are about to construct a OptimalExperimentalDesign java object. First, let R garbage collect
	#to clean up previous objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("OptimalExperimentalDesign.OptimalExperimentalDesign")
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setNandP", as.integer(n), as.integer(p))
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}	
	
	#feed in the data
	for (i in 1 : n){		
		.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
	}
	
	#feed in the inverse var-cov matrix
	if (p <= n){
		for (j in 1 : p){
			.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvXstd[j, , drop = FALSE]) #java indexes from 0...n-1
		}
	}
		
	#now return information as an object (just a list)
	optimal_experimental_design_search = list()
	optimal_experimental_design_search$start = start
	optimal_experimental_design_search$wait = wait
	optimal_experimental_design_search$X = X
	optimal_experimental_design_search$n = n
	optimal_experimental_design_search$p = p
	optimal_experimental_design_search$objective = objective
	optimal_experimental_design_search$java_obj = java_obj
	class(optimal_experimental_design_search) = "optimal_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startOptimalSearch(optimal_experimental_design_search)
	}
	#return the final object
	optimal_experimental_design_search
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param x			The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print optimal_experimental_design_search
#' @export
print.optimal_experimental_design_search = function(x, ...){
	progress = optimalSearchCurrentProgress(x)
	time_elapsed = optimalSearchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the OptimalExperimentalDesign. Did you run \"startOptimalSearch?\"\n")
	} else if (progress == x$max_designs){
		cat("The search completed in", time_elapsed, "seconds.", progress, "vectors have been found.\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%) in ", time_elapsed," seconds.\n", sep = "")
	}
}

#' Prints a summary of a \code{optimal_experimental_design_search} object
#' 
#' @param object		The \code{optimal_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_experimental_design_search
#' @export
summary.optimal_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Plots an order statistic of the object value as a function of number of searches
#' 
#' @param obj			The \code{optimal_experimental_design_search} object whose search history is to be visualized
#' @param order_stat 	The order statistic that you wish to plot. The default is \code{1} for the minimum.
#' @param skip_every	Plot every nth point. This makes the plot generate much more quickly. The default is \code{5}.
#' @param type			The type parameter for plot.
#' @param ... 			Other arguments to be passed to the plot function.
#' @return 				An array of order statistics as a list element
#' 
#' @author 				Adam Kapelner
#' @export
plot_obj_val_order_statistic = function(obj, order_stat = 1, skip_every = 5, type = "o", ...){
	progress = optimalSearchCurrentProgress(obj)
	res = resultsOptimalSearch(ged, max_vectors = 2)	
	vals = res$obj_vals_orig_order
	val_order_stats = array(NA, progress)
	for (d in order_stat : progress){
		if (d %% skip_every == 0){
			val_order_stats[d] = ifelse(order_stat == 1, min(vals[1 : d]), sort(vals[1 : d])[order_stat])
		}		
	}
	plot(1 : progress, val_order_stats, 
			xlab = "Number of Searches", 
			ylab = paste("objective value (", order_stat, ")", sep = ""), 
			type = type, ...)
	invisible(list(val_order_stats = val_order_stats))	
}

#' Starts the parallelized optimal design search. Once begun, this function cannot be run again.
#' 
#' @param obj 		The \code{optimal_experimental_design} object that will be running the search
#' 
#' @author Adam Kapelner
#' @export
startOptimalSearch = function(obj){
	if (.jcall(obj$java_obj, "Z", "began")){
		stop("Search Already begun.")
	}
	.jcall(obj$java_obj, "V", "beginSearch")
}

#' Stops the parallelized optimal design search. Once stopped, it cannot be restarted.
#' 
#' @param obj 		The \code{optimal_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
stopOptimalSearch = function(obj){
	.jcall(obj$java_obj, "V", "stopSearch")
}

# Returns the number of vectors found by the optimal design search
# 
# @param obj 		The \code{optimal_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
optimalSearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}

# Returns the number of vectors found by the optimal design search
# 
# @param obj 		The \code{optimal_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
optimalSearchTimeElapsed = function(obj){
	.jcall(obj$java_obj, "I", "timeElapsedInSeconds")
}

#' Returns the results (thus far) of the optimal design search
#' 
#' @param obj 			The \code{optimal_experimental_design} object that is currently running the search
#' @param max_vectors	The number of design vectors you wish to return. \code{NULL} returns all of them. 
#' 						This is not recommended as returning over 1,000 vectors is time-intensive. The default is 5. 
#' 
#' @author Adam Kapelner
#' @export
resultsOptimalSearch = function(obj, max_vectors = 5){
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	num_iters = .jcall(obj$java_obj, "[I", "getNumIters")
	#these two are in order, so let's order the indicTs by the final objective values
	ordered_indices = order(obj_vals)
	last_index = ifelse(is.null(max_vectors), obj$max_designs, min(max_vectors, obj$max_designs))
	
	ending_indicTs = NULL
	starting_indicTs = NULL
	switches = NULL
	xbarj_diffs = NULL
	pct_vec_same = NULL
	if (max_vectors > 0){
		ending_indicTs = sapply(.jcall(obj$java_obj, "[[I", "getEndingIndicTs", as.integer(ordered_indices[1 : last_index] - 1)), .jevalArray)
	}
	list(
		obj_vals = obj_vals[ordered_indices], 
		num_iters = num_iters[ordered_indices], 
		orig_order = ordered_indices, 
		ending_indicTs = ending_indicTs,
		starting_indicTs = starting_indicTs,
		pct_vec_same = pct_vec_same,
		switches = switches,
		xbarj_diffs = xbarj_diffs
	)
}
