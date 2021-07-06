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
#' @param objective			The method used to aggregate the kernel objective functions together. Default is "added_pct_reduction".
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
#' @param kernel_weights	A vector with positive weights (need not be normalized) where each element represents the weight of 
#' 							each kernel. The default is \code{NULL} for uniform weighting.
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
		objective = "added_pct_reduction",
		kernel_pre_num_designs = 2000,
		kernel_names = NULL,
		Kgrams = NULL,
		maximum_gain_scaling = 1.1,
		kernel_weights = NULL,
		wait = FALSE, 
		start = TRUE,
		max_iters = Inf,
		semigreedy = FALSE, 
		diagnostics = FALSE,
		num_cores = 1,
		seed = NULL){
	
	assertMatrix(X)
	Xstd = standardize_data_matrix(X)
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
	
	n = nrow(Xstd)
	p = ncol(Xstd)
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
						xi = Xstd[i, , drop = FALSE]
						xj = Xstd[j, , drop = FALSE]
						
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
	
	if (!is.null(kernel_weights)){
		assertNumeric(kernel_weights, lower = 0, len = m)
		#normalize here
		kernel_weights = kernel_weights / sum(kernel_weights)
	} else {
		kernel_weights = rep(1 / m, m)
	}
	
	#we are about to construct a GreedyExperimentalDesign java object. First, let R garbage collect
	#to clean up previous GreedyExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!
	
	
	if (objective %in% c("added_pct_reduction")){
		#In this method, we explore each kernel case-by-case before doing a search using all of them
		all_starting_vecs = list()
		all_univariate_kernel_data = list()
		if (is.null(seed)){
			prelim_seed = sample.int(1e5, 1)
		} else {
			prelim_seed = seed
		}
		for (i_k in 1 : m){		
			gd = suppressWarnings(initGreedyExperimentalDesignObject(Xstd, 
							max_designs = kernel_pre_num_designs, Kgram = Kgrams[[i_k]], objective = "kernel",
							diagnostics = TRUE, wait = TRUE, num_cores = nC, seed = prelim_seed)) #same seed should guarantee same starting vectors
			gd_res = resultsGreedySearch(gd, max_vectors = kernel_pre_num_designs, form = "pos_one_min_one")
			objvalsi = array(NA, kernel_pre_num_designs)
			for (i in 1 : kernel_pre_num_designs){
				objvalsi[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
			}
			obj_val_by_iters = gd_res$obj_val_by_iters
			
			#now we want to add
			
			univariate_kernel_data = data.frame(
				kernel = i_k, 
				log10objvalsi = log10(rep(objvalsi, gd_res$num_iters + 1)), 
				log10objvalsf = log10(unlist(gd_res$obj_val_by_iters))
			)
			#we want to eliminate
			univariate_kernel_data = subset(univariate_kernel_data, abs(log10objvalsi - log10objvalsf) > 1e-9)
			univariate_kernel_data$log10_i_over_f = univariate_kernel_data$log10objvalsi - univariate_kernel_data$log10objvalsf
			univariate_kernel_data = univariate_kernel_data[]
			univariate_kernel_data$pct_red_max = 1 - univariate_kernel_data$log10_i_over_f / 
					(max(univariate_kernel_data$log10_i_over_f) * maximum_gain_scaling)
			all_univariate_kernel_data[[i_k]] = univariate_kernel_data
			all_starting_vecs[[i_k]] = gd_res$starting_indicTs
		}		
	}


	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("MultipleKernelGreedyExperimentalDesign.MultipleKernelGreedyExperimentalDesign")
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
	.jcall(java_obj, "V", "setObjective", objective)
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (max_iters < Inf){
		.jcall(java_obj, "V", "setMaxIters", as.integer(max_iters))
	}
#	#feed in the raw data
#	for (i in 1 : n){	
#		.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), Xstd[i, , drop = FALSE]) #java indexes from 0...n-1
#	}
	#feed in the raw kernel data
	for (i_k in 1 : m){
		for (i in 1 : n){	
			.jcall(java_obj, "V", "setSpecificKgramByRow", as.integer(i_k - 1), as.integer(i - 1), Kgrams[[i_k]][i, , drop = FALSE]) #java indexes from 0...n-1
		}	
	}
	
	if (objective %in% c("added_pct_reduction")){
		#feed in data from each kernel in isolation
		for (i_k in 1 : m){
			#.jcall(java_obj, "V", "setObjValLibrary", as.integer(i_k - 1), all_univariate_kernel_data[[i_k]]$log10_i_over_f)
			.jcall(java_obj, "V", "setMaxReductionLogObjVal", as.integer(i_k - 1), max(all_univariate_kernel_data[[i_k]]$log10_i_over_f))
		}
		#set magnification factor
		.jcall(java_obj, "V", "setMaximumGainScaling", maximum_gain_scaling)
		#feed in weights
		.jcall(java_obj, "V", "setKernelWeights", kernel_weights)
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
	ged = list()
	ged$kernel_names = kernel_names
	ged$Kgrams = Kgrams
	ged$m = m
	ged$all_univariate_kernel_data = all_univariate_kernel_data
	ged$kernel_weights = kernel_weights
	ged$all_starting_vecs = all_starting_vecs
	ged$max_designs = max_designs
	ged$kernel_pre_num_designs = kernel_pre_num_designs
	ged$maximum_gain_scaling = maximum_gain_scaling
	ged$semigreedy = semigreedy
	ged$start = start
	ged$wait = wait
	ged$diagnostics = diagnostics
	ged$X = X
	ged$Xstd = Xstd
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

#' Returns the results (thus far) of the greedy design search for multiple kernels
#' 
#' @param obj 			The \code{greedy_multiple_kernel_experimental_design} object that is currently running the search
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
#' 	ged = initGreedyMultipleKernelExperimentalDesignObject(X, 
#' 		max_designs = 1000, num_cores = 3, kernel_names = c("mahalanobis", "gaussian"))
#' 	#wait
#' 	res = resultsMultipleKernelGreedySearch(ged, max_vectors = 2)
#' 	design = res$ending_indicTs[, 1] #ordered already by best-->worst
#'  design
#' 	#how far have we come of the 1000 we set out to do?
#' 	ged
#' 	#we can cut it here
#' 	stopSearch(ged)
#' 	}
#' @export
resultsMultipleKernelGreedySearch = function(obj, max_vectors = 9, form = "one_zero"){
	#get standard information
	greedy_multiple_kernel_experimental_design_search_results = resultsGreedySearch(obj, max_vectors, form)
	last_index = greedy_multiple_kernel_experimental_design_search_results$last_index
	ordered_indices = greedy_multiple_kernel_experimental_design_search_results$ordered_indices
	#now get information specific to this procedure
	if (obj$diagnostics){
		greedy_multiple_kernel_experimental_design_search_results$kernel_obj_vals_by_iter = 
			.jcall(obj$java_obj, "[[[D", "getObjValuesByKernel", as.integer(ordered_indices[1 : last_index] - 1), simplify = TRUE)
	}
	class(greedy_multiple_kernel_experimental_design_search_results) = "greedy_multiple_kernel_experimental_design_search_results"
	#return the final object
	greedy_multiple_kernel_experimental_design_search_results
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