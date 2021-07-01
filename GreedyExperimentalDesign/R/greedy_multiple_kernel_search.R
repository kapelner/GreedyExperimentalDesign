#' Begin A Greedy Pair Multiple Kernel Switching Search
#' 
#' This method creates an object of type greedy_multiple_kernel_experimental_design and will immediately initiate
#' a search through $1_{T}$ space for forced balance designs. For debugging, you can use set the \code{seed}
#' parameter and \code{num_cores = 1} to be assured of deterministic output.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design. We will standardize this matrix by column internally.
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopSearch}} method
#' @param kernel_pre_num_designs	How many designs per kernel to run to explore the space of kernel objective values. Default is 2000.
#' @param kernel_names		An array with the kernels to compute with default parameters. Must have elements in the following set:
#' 							"mahalanobis", "poly_s" where the "s" is a natural number 1 or greater,
#' 							"exponential", "laplacian", "inv_mult_quad", "gaussian". Default is \code{NULL} to
#' 							indicate the kernels are specified manually using the \code{Kgrams} parameter.
#' @param Kgrams			A list of M >= 1 elements where each is a \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. 
#' 							Default is \code{NULL} to indicate this was specified using the convenience parameter 
#' 							\code{kernel_names}.
#' @param maximum_gain_scaling	This controls how much the percentage of possible improvement on a kernel objective function
#' 								should be scaled by. The minimum is 1 which allows for designs that could potentially have >=100%
#' 								improvement over original. We recommend 1.1 which means that a design that was found to be the best
#' 								of the \code{kernel_pre_num_designs} still has 1/1.1 = 9\% room to grow making it highly unlikely that
#' 								any design could be >= 100\%.
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
#' 	ged = initGreedyMultipleKernelExperimentalDesignObject(X, 
#' 		max_designs = 1000, num_cores = 3, kernel_names = c("mahalanobis", "gaussian"))
#' 	#wait
#' 	ged
#' 	}
#' @export
initGreedyMultipleKernelExperimentalDesignObject = function(
		X = NULL, 
		max_designs = 10000, 
		kernel_pre_num_designs = 2000,
		kernel_names = NULL,
		Kgrams = NULL,
		maximum_gain_scaling = 1.1,
		wait = FALSE, 
		start = TRUE,
		max_iters = Inf,
		semigreedy = FALSE, 
		diagnostics = FALSE,
		num_cores = 1,
		seed = NULL){
	
	assertMatrix(X)
	if ((is.null(kernel_names) & is.null(Kgrams)) | (!is.null(kernel_names) & !is.null(Kgrams))){
		stop("You must specify EITHER the kernel_names or the Kgrams argument.")
	}
	assertCount(max_designs, positive = TRUE)
	assertCount(kernel_pre_num_designs, positive = TRUE)
	if (!is.infinite(max_iters)){
		assertCount(max_iters, positive = TRUE)
	}
	assertLogical(wait)
	assertLogical(start)
	assertLogical(semigreedy)
	assertLogical(diagnostics)
	assertCount(num_cores, positive = TRUE)
	assertNumeric(seed, null.ok = TRUE)
	assertNumeric(maximum_gain_scaling, lower = 1)
	
	n = nrow(X)
	p = ncol(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	
	if (!is.null(kernel_names)){
		Kgrams = list()
		for (i_k in 1 : length(kernel_names)){
			if (kernel_names[i_k] == "mahalanobis"){
				Kgrams[[i_k]] = X %*% solve(var(X)) %*% t(X)
			} else {
				Kgrams[[i_k]] = matrix(NA, n, n)
				for (i in 1 : n){
					for (j in 1 : n){
						xi = X[i, , drop = FALSE]
						xj = X[j, , drop = FALSE]
						
						if (str_detect(kernel_names[i_k], "^poly_")){
							s = as.integer(stri_replace_first_regex(kernel_names[i_k], "^poly_(\\d+)", "$1"))
							Kgrams[[i_k]][i, j] = (1 + xi %*% t(xj) / s)^s
						} else if (kernel_names[i_k] == "exponential"){
							Kgrams[[i_k]][i, j] = exp(xi %*% t(xj))
						} else if (kernel_names[i_k] == "laplacian"){
							Kgrams[[i_k]][i, j] = exp(-sqrt(sum((xi - xj)^2)))
						} else if (kernel_names[i_k] == "inv_mult_quad"){
							Kgrams[[i_k]][i, j] = 1 / sqrt(sum((xi - xj)^2) + 1)
						} else if (kernel_names[i_k] == "gaussian"){
							Kgrams[[i_k]][i, j] = exp(-sum((xi - xj)^2))
						} else {
							stop(paste0("Invalid kernel name: ", kernel_names[i_k]))
						}
					}
				}
			}
			
		}
	}
	m = length(Kgrams)
	
	for (i_k in 1 : m){
		assertMatrix(Kgrams[[i_k]], nrows = n, ncols = n)
	}
	
	#we are about to construct a GreedyExperimentalDesign java object. First, let R garbage collect
	#to clean up previous GreedyExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!
	
	#Now we explore each kernel case-by-case
	all_starting_vecs = list()
	all_univariate_kernel_data = data.frame()
	for (i_k in 1 : m){		
		gd = suppressWarnings(initGreedyExperimentalDesignObject(X, 
				max_designs = kernel_pre_num_designs, Kgram = Kgrams[[i_k]], objective = "kernel",
				diagnostics = TRUE, wait = TRUE, num_cores = nC, seed = 1)) #same seed should guarantee same starting vectors
		gd_res = resultsGreedySearch(gd, max_vectors = kernel_pre_num_designs, form = "pos_one_min_one")
		objvalsi = array(NA, Nw)
		for (i in 1 : Nw){
			objvalsi[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
		}
		objvalsf = gd_res$obj_vals_unordered
		
		univariate_kernel_data = data.frame(kernel = i_k, log10objvalsi = log10(objvalsi), log10objvalsf = log10(objvalsf))
		univariate_kernel_data$log10_i_over_f = univariate_kernel_data$log10objvalsi - univariate_kernel_data$log10objvalsf
		univariate_kernel_data$pct_red_max = 1 - univariate_kernel_data$log10_i_over_f / 
													(max(univariate_kernel_data$log10_i_over_f) * maximum_gain_scaling)
		all_univariate_kernel_data = rbind(all_univariate_kernel_data, univariate_kernel_data)
		all_starting_vecs[[i_k]] = gd_res$starting_indicTs
	}

	
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
	.jcall(java_obj, "V", "setN", as.integer(n))
	.jcall(java_obj, "V", "setP", as.integer(p))
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (max_iters < Inf){
		.jcall(java_obj, "V", "setMaxIters", as.integer(max_iters))
	}
	
	
	
	
	
	
	
	#feed in the gram matrix if applicable
#	if (!is.null(Kgram)){
#		setGramMatrix(java_obj, Kgram)
#	} else {
#		#feed in the raw data
#		for (i in 1 : n){	
#			if (objective == "abs_sum_diff"){
#				.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
#			} else {
#				.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
#			}
#		}
#		
#		#feed in the inverse var-cov matrix
#		if (objective == "mahal_dist"){
#			if (p < n){
#				for (j in 1 : p){
#					.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
#				}
#			}
#		}
#	}
	
	#do we want diagnostics? Set it...
	if (diagnostics){
		.jcall(java_obj, "V", "setDiagnostics")
	}
	
	#is it semigreedy? Set it...
	if (semigreedy){
		.jcall(java_obj, "V", "setSemigreedy")
	}
	
	#now return information as an object (just a list)
	ged = list()
	ged$kernel_names = kernel_names
	ged$Kgrams = Kgrams
	ged$m = m
	ged$all_univariate_kernel_data = all_univariate_kernel_data
	ged$all_starting_vecs = all_starting_vecs
	ged$max_designs = max_designs
	ged$kernel_pre_num_designs = kernel_pre_num_designs
	ged$maximum_gain_scaling = maximum_gain_scaling
	ged$semigreedy = semigreedy
	ged$start = start
	ged$wait = wait
	ged$diagnostics = diagnostics
	ged$X = X
	ged$n = n
	ged$p = p
	ged$java_obj = java_obj
	class(ged) = "greedy_multiple_kernel_experimental_design"
	#if the user wants to run it immediately...
	if (start){
		startSearch(ged)
	}
	#return the final object
	ged
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
		#we should make switches into a list now
		xbarj_diffs = .jcall(obj$java_obj, "[[[D", "getXbarjDiffs", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		obj_val_by_iters = .jcall(obj$java_obj, "[[D", "getObjValByIter", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
		
		pct_vec_same = colSums(starting_indicTs == ending_indicTs) / length(starting_indicTs[, 1]) * 100
	}
	greedy_experimental_design_search_results = list(
			obj_vals = obj_vals[ordered_indices], 
			obj_vals_unordered = obj_vals,
			num_iters = num_iters[ordered_indices], 
			orig_order = ordered_indices, 
			ending_indicTs = ending_indicTs,
			starting_indicTs = starting_indicTs,
			obj_val_by_iters = obj_val_by_iters,
			pct_vec_same = pct_vec_same,
			switches = switches,
			xbarj_diffs = xbarj_diffs
	)
	class(greedy_experimental_design_search_results) = "greedy_experimental_design_search_results"
	#return the final object
	greedy_experimental_design_search_results
}


#' Prints a summary of a \code{greedy_multiple_kernel_experimental_design} object
#' 
#' @param x			The \code{greedy_multiple_kernel_experimental_design} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print greedy_multiple_kernel_experimental_design
#' @export
print.greedy_multiple_kernel_experimental_design = function(x, ...){
	print.greedy_experimental_design_search(x, ...)
}

#' Prints a summary of a \code{greedy_multiple_kernel_experimental_design} object
#' 
#' @param object		The \code{greedy_multiple_kernel_experimental_design} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary greedy_multiple_kernel_experimental_design
#' @export
summary.greedy_multiple_kernel_experimental_design = function(object, ...){
	print(object, ...)
}

#' Plots a summary of a \code{greedy_multiple_kernel_experimental_design} object
#' 
#' @param x			The \code{greedy_multiple_kernel_experimental_design} object to be summarized in the plot
#' @param ...		Other parameters to pass to the default plot function
#' @return			An array of order statistics from \link{plot_obj_val_order_statistic} as a list element
#' 
#' @author 			Adam Kapelner
#' @method plot greedy_experimental_design_search
#' @export
plot.greedy_multiple_kernel_experimental_design = function(x, ...){
	plot.greedy_experimental_design_search(x, ...)
}