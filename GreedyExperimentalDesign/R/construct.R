VERSION = "1.0"

#' This method creates an object of type greedy_experimental_design and will immediately initiate
#' a search through $1_{T}$ space.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopGreedySearch}} method 
#' @param objective			The objective function to use when greedily searching design space. This is a string
#' 							"\code{abs_sum_diff}" (default) or "\code{mahal_dist}."
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @param semigreedy		Should we use a fully greedy approach or the quicker semi-greedy approach? The default is
#' 							\code{FALSE} corresponding to the fully greedy approach.
#' @param max_iters			Should we impose a maximum number of greedy switches? The default is \code{Inf} which a flag 
#' 							for ``no limit.''
#' @param diagnostics		Returns diagnostic information about the iterations including (a) the initial starting
#' 							vectors, the switches at every iteration and information about the objective function
#' 							at every iteration (default is \code{FALSE} due to speed concerns).
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @return					An object of type \code{greedy_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initGreedyExperimentalDesignObject = function(X, 
		max_designs = 10000, 
		objective = "abs_sum_diff", 
		wait = FALSE, 
		start = TRUE,
		max_iters = Inf,
		semigreedy = FALSE, 
		diagnostics = FALSE,
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
	
	#we are about to construct a GreedyExperimentalDesign java object. First, let R garbage collect
	#to clean up previous GreedyExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("GreedyExperimentalDesign.GreedyExperimentalDesign")
	.jcall(java_obj, "V", "setMaxDesigns", as.integer(max_designs))
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
	.jcall(java_obj, "V", "setNandP", as.integer(n), as.integer(p))
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (max_iters <= 0){stop("max_iters must be positive")}
	if (max_iters < Inf){
		.jcall(java_obj, "V", "setMaxIters", as.integer(max_iters))
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
	
	#do we want diagnostics? Set it...
	if (diagnostics){
		.jcall(java_obj, "V", "setDiagnostics")
	}
	
	#is it semigreedy? Set it...
	if (semigreedy){
		.jcall(java_obj, "V", "setSemigreedy")
	}
		
	#now return information as an object (just a list)
	greedy_experimental_design_search = list()
	greedy_experimental_design_search$max_designs = max_designs
	greedy_experimental_design_search$semigreedy = semigreedy
	greedy_experimental_design_search$X = X
	greedy_experimental_design_search$n = n
	greedy_experimental_design_search$p = p
	greedy_experimental_design_search$objective = objective
	greedy_experimental_design_search$java_obj = java_obj
	class(greedy_experimental_design_search) = "greedy_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startGreedySearch(greedy_experimental_design_search)
	}
	#return the final object
	greedy_experimental_design_search
}

#' Prints a summary of a \code{greedy_experimental_design_search} object
#' 
#' @param x			The \code{greedy_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print greedy_experimental_design_search
#' @export
print.greedy_experimental_design_search = function(x, ...){
	progress = greedySearchCurrentProgress(x)
	time_elapsed = greedySearchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the GreedyExperimentalDesign. Did you run \"startGreedySearch?\"\n")
	} else if (progress == x$max_designs){
		cat("The search completed in", time_elapsed, "seconds.", progress, "vectors have been found.\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%) in ", time_elapsed," seconds.\n", sep = "")
	}
}

#' Prints a summary of a \code{greedy_experimental_design_search} object
#' 
#' @param object		The \code{greedy_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary greedy_experimental_design_search
#' @export
summary.greedy_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Plots a summary of a \code{greedy_experimental_design_search} object
#' 
#' @param x			The \code{greedy_experimental_design_search} object to be summarized in the plot
#' @param ...		Other parameters to pass to the default plot function
#' @return			An array of order statistics from \link{\code{plot_obj_val_order_statistic}} as a list element
#' 
#' @author 			Adam Kapelner
#' @method plot greedy_experimental_design_search
#' @export
plot.greedy_experimental_design_search = function(x, ...){
	par(mfrow = c(1, 2))
	
	progress = greedySearchCurrentProgress(x)
	res = resultsGreedySearch(ged, max_vectors = 2)
	hist(res$obj_vals_orig_order, br = progress / 10, xlab = "objective value", ylab = NULL, main = paste("After", progress, "searches"))
#	hist(res$num_iters, br = progress / 10, xlab = "# of search iterations", ylab = NULL, main = "")
	
	#now do the plot of number of searches needed
	plot_obj_val_order_statistic(x)
}

#' Plots an order statistic of the object value as a function of number of searches
#' 
#' @param obj			The \code{greedy_experimental_design_search} object whose search history is to be visualized
#' @param order_stat 	The order statistic that you wish to plot. The default is \code{1} for the minimum.
#' @param skip_every	Plot every nth point. This makes the plot generate much more quickly. The default is \code{5}.
#' @param type			The type parameter for plot.
#' @param ... 			Other arguments to be passed to the plot function.
#' @return 				An array of order statistics as a list element
#' 
#' @author 				Adam Kapelner
#' @export
plot_obj_val_order_statistic = function(obj, order_stat = 1, skip_every = 5, type = "o", ...){
	progress = greedySearchCurrentProgress(obj)
	res = resultsGreedySearch(ged, max_vectors = 2)	
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

#' Starts the parallelized greedy design search. Once begun, this function cannot be run again.
#' 
#' @param obj 		The \code{greedy_experimental_design} object that will be running the search
#' 
#' @author Adam Kapelner
#' @export
startGreedySearch = function(obj){
	if (.jcall(obj$java_obj, "Z", "began")){
		stop("Search Already begun.")
	}
	.jcall(obj$java_obj, "V", "beginSearch")
}

#' Stops the parallelized greedy design search. Once stopped, it cannot be restarted.
#' 
#' @param obj 		The \code{greedy_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
stopGreedySearch = function(obj){
	.jcall(obj$java_obj, "V", "stopSearch")
}

# Returns the number of vectors found by the greedy design search
# 
# @param obj 		The \code{greedy_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
greedySearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}

# Returns the number of vectors found by the greedy design search
# 
# @param obj 		The \code{greedy_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
greedySearchTimeElapsed = function(obj){
	.jcall(obj$java_obj, "I", "timeElapsedInSeconds")
}

#' Returns the results (thus far) of the greedy design search
#' 
#' @param obj 			The \code{greedy_experimental_design} object that is currently running the search
#' @param max_vectors	The number of design vectors you wish to return. \code{NULL} returns all of them. 
#' 						This is not recommended as returning over 1,000 vectors is time-intensive. The default is 5. 
#' 
#' @author Adam Kapelner
#' @export
resultsGreedySearch = function(obj, max_vectors = 5){
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	num_iters = .jcall(obj$java_obj, "[I", "getNumIters")
	#these two are in order, so let's order the indicTs by the final objective values
	ordered_indices = order(obj_vals)
	last_index = ifelse(is.null(max_vectors), length(ordered_indices), max_vectors)
	
	indicTs = NULL
	if (max_vectors > 0){
		indicTs = sapply(.jcall(obj$java_obj, "[[I", "getEndingIndicTs", as.integer(ordered_indices[1 : last_index] - 1)), .jevalArray)
	}
	list(obj_vals = obj_vals[ordered_indices], num_iters = num_iters[ordered_indices], obj_vals_orig_order = obj_vals, indicTs = indicTs)
}

# PRIVATE
compute_objectives = function(X, indic_T, inv_cov_X = NULL){
	#saves computation sometimes to pass it in
	if (is.null(inv_cov_X)){
		inv_cov_X = solve(var(X))
	}
	
	X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_C = X[indic_T == 0, , drop = FALSE] #coerce as matrix in order to matrix multiply later
	X_T_bar = colMeans(X_T)
	X_C_bar = colMeans(X_C)
	X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar)
	
	abs_obj = sum(abs(X_T_bar_minus_X_C_bar))
	#do the main calculation (eq 6) and return it as a scalar
	
	mahal_obj = as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)	
	
	list(abs_obj = abs_obj, mahal_obj = mahal_obj)
}


#' Generates a design matrix with standardized predictors. Useful for debugging.
#' 
#' @param n					Number of rows in the design matrix 
#' @param p 				Number of columns in the design matrix
#' @param covariate_dist	The function to use to draw the covariate realizations (assumed to be iid).
#' 							This defaults to \code{rnorm} for $N(0,1)$ draws.
#' @param ...				Optional arguments to be passed to the \code{covariate_dist} function.
#' @return 					THe design matrix
#' 
#' @author Adam Kapelner
#' @export
generate_stdzied_design_matrix = function(n = 50, p = 1, covariate_gen = rnorm, ...){
	X = matrix(covariate_gen(n * p, ...), nrow = n, ncol = p)
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}

# PRIVATE: Creates a random binary vector which codes an experimental design
# 
# @param n			How many subjects does the experiment have 
# @param p_w 		What is the proportion of treatments?
# @return 			An array of length $n$ with $n \times p_w$ number of 1's and the rest 0's in a random order
# 
# @author Adam Kapelner
#create_random_dummy_vec = function(n, p_w = 0.5){
#	indic_T_dummy_permutation_vector = c(rep(1, n * p_w), rep(0, n * p_w)) #there are n_T 1's followed by n_C 0's dictated by p_w
#	sample(indic_T_dummy_permutation_vector)
#}


