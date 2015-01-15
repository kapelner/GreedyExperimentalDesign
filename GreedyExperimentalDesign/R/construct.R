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
#' @param num_cores 		The number of CPU cores you wish to use during the search
#' @return					An object of type \code{greedy_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initGreedyExperimentalDesignObject = function(X, max_designs = 10000, objective = "abs_sum_diff", num_cores = 1){
	#get dimensions immediately
	n = nrow(X)
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
	
	#now feed into Java some starting points for the search since it's easier to do in R
	for (d in 1 : max_designs){
		.jcall(java_obj, "V", "setDesignStartingPoint", as.integer(d - 1), as.integer(create_random_dummy_vec(n))) #java indexes from 0...n-1
	}
		
	#now return information as an object (just a list)
	greedy_experimental_design_search = list()
	greedy_experimental_design_search$max_designs = max_designs
	greedy_experimental_design_search$X = X
	greedy_experimental_design_search$n = n
	greedy_experimental_design_search$p = p
	greedy_experimental_design_search$objective = objective
	greedy_experimental_design_search$java_obj = java_obj
	class(greedy_experimental_design_search) = "greedy_experimental_design_search"
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
	if (progress == 0){
		cat("No progress on the GreedyExperimentalDesign. Did you run \"startGreedySearch?\"\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%).\n", sep = "")
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
#' 
#' @author 			Adam Kapelner
#' @method plot greedy_experimental_design_search
#' @export
plot.greedy_experimental_design_search = function(x, ...){
	par(mfrow = c(1, 2))
	
	progress = greedySearchCurrentProgress(x)
	res = resultsGreedySearch(ged, max_vectors = 2)
	hist(res$obj_vals_orig_order, br = progress / 10, xlab = "objective value", ylab = NULL, main = paste("After", progress, "searches"))
	
	obj_vals_rand_order = res$obj_vals_orig_order
	obj_vals_rand_order_mins = array(NA, progress)
	for (d in 1 : progress){
		obj_vals_rand_order_mins[d] = min(obj_vals_rand_order[1 : d])
	}
	plot(1 : progress, obj_vals_rand_order_mins, xlab = "Number of Searches", ylab = "Best objective value", ...)
	invisible(list(obj_vals_rand_order_mins = obj_vals_rand_order_mins))
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

#' Returns the number of vectors found by the greedy design search
#' 
#' @param obj 		The \code{greedy_experimental_design} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @export
greedySearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}

#' Returns the results (thus far) of the greedy design search
#' 
#' @param obj 			The \code{greedy_experimental_design} object that is currently running the search
#' @param max_vectors	The number of design vectors you wish to return. \code{NULL} returns all of them 
#' 						(default). This is not recommended as returning over 1,000 vectors is time-intensive.
#' 
#' @author Adam Kapelner
#' @export
resultsGreedySearch = function(obj, max_vectors = NULL){
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	#these two are in order, so let's order the indicTs by the final objective values
	ordered_indices = order(obj_vals)
	last_index = ifelse(is.null(max_vectors), length(ordered_indices), max_vectors)
	
	
	indicTs = t(sapply(.jcall(obj$java_obj, "[[I", "getEndingIndicTs", as.integer(ordered_indices[1 : last_index])), .jevalArray))
	list(obj_vals = obj_vals[ordered_indices], obj_vals_orig_order = obj_vals, indicTs = indicTs)
}

# PRIVATE: Creates a random binary vector which codes an experimental design
# 
# @param n			How many subjects does the experiment have 
# @param p_w 		What is the proportion of treatments?
# @return 			An array of length $n$ with $n \times p_w$ number of 1's and the rest 0's in a random order
# 
# @author Adam Kapelner
create_random_dummy_vec = function(n, p_w = 0.5){
	indic_T_dummy_permutation_vector = c(rep(1, n * p_w), rep(0, n * p_w)) #there are n_T 1's followed by n_C 0's dictated by p_w
	sample(indic_T_dummy_permutation_vector)
}


