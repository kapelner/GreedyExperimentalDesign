#' Begin A Greedy Pair Switching Search
#' 
#' This method creates an object of type greedy_experimental_design and will immediately initiate
#' a search through $1_{T}$ space for forced balance designs. For debugging, you can use set the \code{seed}
#' parameter and \code{num_cores = 1} to be assured of deterministic output.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design. This parameter must be specified unless you
#' 							choose objective type \code{"kernel"} in which case, the \code{Kgram} parameter must
#' 							be specified.
#' @param nT				The number of treatments to assign. Default is \code{NULL} which is for forced balance allocation
#' 							i.e. nT = nC = n / 2 where n is the number of rows in X (or Kgram if X is unspecified).
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopSearch}} method 
#' @param objective			The objective function to use when searching design space. This is a string
#' 							with valid values "\code{mahal_dist}" (the default), "\code{abs_sum_diff}" or "\code{kernel}".
#' @param indicies_pairs	A matrix of size $n/2$ times 2 whose rows are indicies pairs. The values of the entire matrix 
#' 							must enumerate all indicies $1, ..., n$. The default is \code{NULL} meaning to use all possible pairs.
#' @param Kgram				If the \code{objective = kernel}, this argument is required to be an \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param wait				Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 							deafult is \code{FALSE}.
#' @param start				Should we start searching immediately (default is \code{TRUE}).
#' @param semigreedy		Should we use a fully greedy approach or the quicker semi-greedy approach? The default is
#' 							\code{FALSE} corresponding to the fully greedy approach.
#' @param max_iters			Should we impose a maximum number of greedy switches? The default is \code{Inf} which a flag 
#' 							for ``no limit.''
#' @param diagnostics		Returns diagnostic information about the iterations including (a) the initial starting
#' 							vectors, (b) the switches at every iteration and (c) information about the objective function
#' 							at every iteration (default is \code{FALSE} to decrease the algorithm's run time).
#' @param num_cores 		The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @param seed				The set to set for deterministic output. This should only be set if \code{num_cores = 1} otherwise
#' 							the output will not be deterministic. Default is \code{NULL} for no seed set.
#' @return					An object of type \code{greedy_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#' 	library(MASS)
#' 	data(Boston)
#'  #pretend the Boston data was an experiment setting 
#' 	#first pull out the covariates
#'  X = Boston[, 1 : 13] 
#'  #begin the greedy design search
#' 	ged = initGreedyExperimentalDesignObject(X, 
#' 		max_designs = 1000, num_cores = 3, objective = "abs_sum_diff")
#' 	#wait
#' 	ged
#' 	}
#' @export
initGreedyExperimentalDesignObject = function(
		X = NULL, 
		nT = NULL,
		max_designs = 10000, 
		objective = "mahal_dist", 
		indicies_pairs = NULL,
		Kgram = NULL,
		wait = FALSE, 
		start = TRUE,
		max_iters = Inf,
		semigreedy = FALSE, 
		diagnostics = FALSE,
		num_cores = 1,
		seed = NULL){
	
	if (!is.null(Kgram)){
		n = nrow(Kgram)
		p = NA
	} else {
		n = nrow(X)
		p = ncol(X)
	}
	if (is.null(nT)){
		nT = n / 2
	}
	verify_objective_function(objective, Kgram, n)
	
	checkCount(nT, positive = TRUE, null.ok = FALSE)
	
	if (nT != n / 2 & !is.null(indicies_pairs)){
		stop("indicies_pairs cannot be specified if nT is not half of n")
	}
	
	if (!is.null(indicies_pairs)){
		if (!(all.equal(sort(c(indicies_pairs)), 1 : n))){
			stop("indicies_pairs must cover all indicies 1, 2, ..., n once each.")
		}
	}
	if (objective == "abs_sum_diff"){
		#standardize it -- much faster here
		Xstd = standardize_data_matrix(X)
	}
	if (objective == "mahal_dist"){
		if (p < n){
			SinvX = solve(var(X))
		}
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
	if (!is.null(seed)){
		.jcall(java_obj, "V", "setSeed", as.integer(seed))
		if (num_cores != 1){
			warning("Setting the seed with multiple cores does not guarantee deterministic output.")
		}		
	}
	if (objective != "kernel"){
		p = ncol(X)
		.jcall(java_obj, "V", "setP", as.integer(p))
	}
	.jcall(java_obj, "V", "setN", as.integer(n))
	.jcall(java_obj, "V", "setNumTreatments", as.integer(nT))
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (max_iters <= 0){stop("max_iters must be positive")}
	if (max_iters < Inf){
		.jcall(java_obj, "V", "setMaxIters", as.integer(max_iters))
	}
	
	if (!is.null(indicies_pairs)){
		for (i in 1 : (n / 2)){	
			.jcall(java_obj, "V", "setLegalPair", as.integer(indicies_pairs[i, ] - 1), as.integer(i - 1)) #java indexes from 0...n-1
		}
	}
	
	#feed in the gram matrix if applicable
	if (!is.null(Kgram)){
		setGramMatrix(java_obj, Kgram)
	} else {
		#feed in the raw data
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
	greedy_experimental_design_search$start = start
	greedy_experimental_design_search$wait = wait
	greedy_experimental_design_search$diagnostics = diagnostics
	greedy_experimental_design_search$X = X
	greedy_experimental_design_search$n = n
	greedy_experimental_design_search$p = p
	greedy_experimental_design_search$objective = objective
	greedy_experimental_design_search$java_obj = java_obj
	class(greedy_experimental_design_search) = "greedy_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(greedy_experimental_design_search)
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
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the GreedyExperimentalDesign. Did you run \"startSearch?\"\n")
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

#' Plots a summary of a greedy search object object
#' 
#' @param x			The greedy search object object to be summarized in the plot
#' @param ...		Other parameters to pass to the default plot function
#' @return			An array of order statistics from \link{plot_obj_val_order_statistic} as a list element
#' 
#' @author 			Adam Kapelner
#' @method plot greedy_experimental_design_search
#' @export
plot.greedy_experimental_design_search = function(x, ...){
	par(mfrow = c(1, 2))
	
	progress = greedySearchCurrentProgress(x)
	res = resultsGreedySearch(x, max_vectors = 2)
	hist(res$obj_vals_orig_order, br = progress / 10, xlab = "objective value", ylab = NULL, main = paste("After", progress, "searches"))
#	hist(res$num_iters, br = progress / 10, xlab = "# of search iterations", ylab = NULL, main = "")
	
	#now do the plot of number of searches needed
	plot_obj_val_order_statistic(x)
}

#' Plots the objective value by iteration
#' 
#' @param res 		Results from a greedy search object
#' @param runs 		A vector of run indices you would like to see plotted (default is to plot the first up to 9)
#' 
#' @author 			Adam Kapelner
#' @export
plot_obj_val_by_iter = function(res, runs = NULL){
	if (is.null(res$obj_val_by_iters)){
		stop("You need to set diagnostics = TRUE on the search object.")
	}
	
	if (is.null(runs)){
		runs = 1 : min(length(res$obj_val_by_iters), 9)
	}
	num_to_plot = length(runs)
	
	par(mfrow = c(ceiling(sqrt(num_to_plot)), ceiling(sqrt(num_to_plot))))
	for (run in runs){
		obj_vals = res$obj_val_by_iters[[run]]
		plot(1 : length(obj_vals), obj_vals, main = paste("Run #", run), xlab = "iteration", ylab = "objective value", type = "o")
	}
} 

#' Plots an order statistic of the object value as a function of number of searches
#' 
#' @param obj			The greedy search object object whose search history is to be visualized
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
	res = resultsGreedySearch(obj, max_vectors = 0)	#don't need any vectors => set it to 0 for speed concerns
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

# Returns the number of vectors found by the greedy design search
# 
# @param obj 		The \code{greedy_experimental_design} object that is currently running the search
# 
# @author Adam Kapelner
greedySearchCurrentProgress = function(obj){
	.jcall(obj$java_obj, "I", "progress")
}


#' Returns the results (thus far) of the greedy design search
#' 
#' @param obj 			The \code{greedy_experimental_design} object that is currently running the search
#' @param max_vectors	The number of design vectors you wish to return. \code{NULL} returns all of them. 
#' 						This is not recommended as returning over 1,000 vectors is time-intensive. The default is 9. 
#' @param form			Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#' 	library(MASS)
#' 	data(Boston)
#'  #pretend the Boston data was an experiment setting 
#' 	#first pull out the covariates
#'  X = Boston[, 1 : 13]
#'  #begin the greedy design search
#' 	ged = initGreedyExperimentalDesignObject(X, 
#' 		max_designs = 1000, num_cores = 2, objective = "abs_sum_diff")
#' 	#wait
#' 	res = resultsGreedySearch(ged, max_vectors = 2)
#' 	design = res$ending_indicTs[, 1] #ordered already by best-->worst
#'  design
#'  #what is the balance on this vector?
#' 	res$obj_vals[1]
#' 	#compute balance explicitly in R to double check
#' 	compute_objective_val(X, design) #same as above
#' 	#how far have we come?
#' 	ged
#' 	#we can cut it here
#' 	stopSearch(ged)
#' 	}
#' @export
resultsGreedySearch = function(obj, max_vectors = 9, form = "one_zero"){
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	num_iters = .jcall(obj$java_obj, "[I", "getNumIters")
	#these two are in order, so let's order the indicTs by the final objective values
	ordered_indices = order(obj_vals)
	last_index = ifelse(is.null(max_vectors), obj$max_designs, min(max_vectors + 1, obj$max_designs))
	
	ending_indicTs = NULL
	starting_indicTs = NULL
	switches = NULL
	xbarj_diffs = NULL
	obj_val_by_iters = NULL
	pct_vec_same = NULL
	ending_indicTs = .jcall(obj$java_obj, "[[I", "getEndingIndicTs", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
	if (form == "pos_one_min_one"){
		ending_indicTs = (ending_indicTs - 0.5) * 2
	}
	if (obj$diagnostics){
		starting_indicTs = .jcall(obj$java_obj, "[[I", "getStartingIndicTs", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		if (form == "pos_one_min_one"){
			starting_indicTs = (starting_indicTs - 0.5) * 2
		}
		switches = .jcall(obj$java_obj, "[[[I", "getSwitchedPairs", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		xbarj_diffs = .jcall(obj$java_obj, "[[[D", "getXbarjDiffs", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		obj_val_by_iters = .jcall(obj$java_obj, "[[D", "getObjValByIter", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		
		pct_vec_same = colSums(starting_indicTs == ending_indicTs) / length(starting_indicTs[, 1]) * 100
	}
	greedy_experimental_design_search_results = list(
		ordered_indices = ordered_indices,
		obj_vals = obj_vals[ordered_indices], 
		obj_vals_unordered = obj_vals,
		num_iters = num_iters[ordered_indices], 
		orig_order = ordered_indices, 
		ending_indicTs = ending_indicTs,
		starting_indicTs = starting_indicTs,
		obj_val_by_iters = obj_val_by_iters,
		pct_vec_same = pct_vec_same,
		switches = switches,
		xbarj_diffs = xbarj_diffs,
		last_index = last_index
	)
	class(greedy_experimental_design_search_results) = "greedy_experimental_design_search_results"
	#return the final object
	greedy_experimental_design_search_results
}

#' Compute Optimal Number of Treatments/Controls
#' 
#' Given a total budget and asymmetric treatment and control costs, calculate the
#' number of treatments and controls that optimize the variance of the estimator. 
#' The number of treatments is rounded up by default.
#' 
#' @param c_treatment		The cost of a treatment assignment. Default is \code{NULL} for symmetric costs.
#' @param c_control			The cost of a control assignment. Default is \code{NULL} for symmetric costs.
#' @param c_total_max		The total cost constraint of any allocation. Either this or \code{n} must be specified. Default is \code{NULL}.
#' @param n					The total cost constraint as specified by the total number of subjects. Either this or \code{c_total} must be 
#' 							specified. Default is \code{NULL}.	
#' 
#' @return					A list with three keys: n, nT, nC plus specified arguments
#' 
#' @author Adam Kapelner
#' @examples
#'  \dontrun{
#' 	optimize_asymmetric_treatment_assignment(n = 100) 
#'  #nT = nC = 50
#' 	optimize_asymmetric_treatment_assignment(n = 100, c_treatment = 2, c_control = 1) 
#'  #nT = 66, nC = 34
#' 	optimize_asymmetric_treatment_assignment(c_total_max = 50, c_treatment = 2, c_control = 1) 
#' 	}
#' @export
optimize_asymmetric_treatment_assignment = function(
		c_treatment = NULL,
		c_control = NULL,
		c_total_max = NULL,
		n = NULL,
		round_up_nT = FALSE
	){
	if ((is.null(c_total_max) & is.null(n)) | (!is.null(c_total_max) & !is.null(n))){
		stop("n xor c_total must be specified.")
	}
	if (!is.null(n)){
		if (is.null(c_treatment) & is.null(c_control)){
			nT = n / 2
			list(n = n, nT = nT, nC = n - nT)
		} else {
			checkNumeric(c_treatment, lower = .Machine$double.xmin, finite = TRUE)
			checkNumeric(c_control,   lower = .Machine$double.xmin, finite = TRUE)
			c_const = c_treatment / c_control
			nT = if (c_treatment > c_control){
				 	ceiling(n / (1 + sqrt(c_const)))
		         } else {
				 	floor(n / (1 + sqrt(c_const)))
		         }
		 	nC = n - nT
		 	c_total = nT * c_treatment + nC * c_control
		 	list(n = n, nT = nT, nC = nC, c_treatment = c_treatment, c_control = c_control, c_total = c_total)
		}
	} else {
		checkNumeric(c_treatment, lower = .Machine$double.xmin, finite = TRUE)
		checkNumeric(c_control,   lower = .Machine$double.xmin, finite = TRUE)
		
		K_const = c_total_max / c_control
		c_const = c_treatment / c_control
		nT = K_const / (c_const + sqrt(c_const))
		nC = K_const - c_const * nT
		nT = floor(nT)
		nC = if (c_const * nT + ceiling(nC) <= K_const){
					ceiling(nC)
				} else {
					floor(nC)
				}
		c_total = nT * c_treatment + nC * c_control
		list(n = nT + nC, nT = nT, nC = nC, c_total = c_total, c_total_max = c_total_max, c_treatment = c_treatment, c_control = c_control)
	}
}


