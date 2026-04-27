#' Begin A Greedy Multiple Kernel Design Search
#' 
#' This method creates an object of type greedy_multiple_kernel_experimental_design and will immediately initiate
#' a search through allocation space for forced balance designs. For debugging, you can use set the \code{seed}
#' parameter and \code{num_cores = 1} to be assured of deterministic output.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design. This parameter must be specified unless you
#' 							choose objective type \code{"kernel"} in which case, the \code{Kgrams} parameter must
#' 							be specified.
#' @param kernel_names		A vector of M >= 1 strings indicating the kernels to be used. Valid values are 
#' 							"\code{mahalanobis}", "\code{gaussian}", "\code{laplacian}", "\code{exponential}", 
#' 							"\code{inv_mult_quad}", "\code{poly_2}", "\code{poly_3}". Or a special value "\code{manual}" to 
#' 							indicate the kernels are specified manually using the \code{Kgrams} parameter.
#' @param Kgrams			A list of M >= 1 elements where each is a \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param kernel_weights	A vector of M weights for each kernel which should sum to 1. Default is \code{NULL} which 
#' 							specifies equal weighting.
#' @param maximum_gain_scaling Should we scale the kernels to have the same maximum gain? Default is \code{TRUE}.
#' @param nT				The number of treatments to assign. Default is \code{NULL} which is for forced balance allocation
#' 							i.e. nT = nC = n / 2 where n is the number of rows in X (or Kgram if X is unspecified).
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopSearch}} method. 
#' @param kernel_pre_num_designs The number of initial designs to search through before starting the greedy search for each kernel.
#' 							Default is \code{1000}.
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
#' @param verbose			Should the algorithm emit progress output? Default is \code{TRUE}.
#' @return					An object of type \code{greedy_multiple_kernel_experimental_design} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
initGreedyMultipleKernelExperimentalDesignObject = function(
		X = NULL, 
		kernel_names = NULL,
		Kgrams = NULL,
		kernel_weights = NULL,
		maximum_gain_scaling = TRUE,
		nT = NULL,
		max_designs = 10000, 
		kernel_pre_num_designs = 1000,
		wait = FALSE, 
		start = TRUE,
		max_iters = Inf,
		semigreedy = FALSE, 
		diagnostics = FALSE,
		num_cores = 1,
		seed = NULL,
		verbose = TRUE){
	
	if ((is.null(kernel_names) & is.null(Kgrams)) | (!is.null(kernel_names) & !is.null(Kgrams))){
		stop("You must specify EITHER the kernel_names or the Kgrams argument.")
	}
	
	if (!is.null(Kgrams)){
		n = nrow(Kgrams[[1]])
		p = NA
	} else {
		n = nrow(X)
		p = ncol(X)
	}
	if (is.null(nT)){
		nT = n / 2
	}
	
	checkCount(nT, positive = TRUE, null.ok = FALSE)
	assertLogical(verbose)
	assertLogical(maximum_gain_scaling)
	
	if (is.null(Kgrams)){
		Xstd = standardize_data_matrix(X)
		SinvX = solve(stats::var(X))
		Kgrams = list()
		for (i_k in 1 : length(kernel_names)){
			if (kernel_names[i_k] == "mahalanobis"){
				Kgrams[[i_k]] = X %*% SinvX %*% t(X)
			} else {
				if (grepl("poly", kernel_names[i_k])){
					s = as.numeric(sub("poly_", "", kernel_names[i_k]))
					Kgrams[[i_k]] = compute_kernel_matrix_gpu(Xstd, kernel = "poly", poly_s = s)
				} else {
					Kgrams[[i_k]] = compute_kernel_matrix_gpu(Xstd, kernel = kernel_names[i_k])
				}
			}
		}
	}
	
	m = length(Kgrams)
	for (i_k in 1 : m){
		assertMatrix(Kgrams[[i_k]], nrows = n, ncols = n)
	}
	
	if (is.null(kernel_weights)){
		kernel_weights = rep(1 / m, m)
	}
	assertNumeric(kernel_weights, len = m, lower = 0, upper = 1)
	if (abs(sum(kernel_weights) - 1) > 1e-8){
		stop("kernel_weights must sum to 1.")
	}

	#we are about to construct a MultipleKernelGreedyExperimentalDesign java object. First, let R garbage collect
	#to clean up previous GreedyExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("MultipleKernelGreedyExperimentalDesign.MultipleKernelGreedyExperimentalDesign")
	set_verbose_if_available(java_obj, verbose)
	
	.jcall(java_obj, "V", "setMaxDesigns", as.integer(max_designs))
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))	
	if (!is.null(seed)){
		.jcall(java_obj, "V", "setSeed", as.integer(seed))
	}
	.jcall(java_obj, "V", "setN", as.integer(n))
	.jcall(java_obj, "V", "setNumTreatments", as.integer(nT))
	.jcall(java_obj, "V", "setM", as.integer(m))
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	if (max_iters <= 0){stop("max_iters must be positive")}
	if (max_iters < Inf){
		.jcall(java_obj, "V", "setMaxIters", as.integer(max_iters))
	}
	
	#find r many starting points
	#first find optimal vectors for each individual kernel
	cat("Finding initial starting points for each kernel individually.\n")
	kernel_obj_values = matrix(NA, nrow = max_designs, ncol = m)
	initial_starting_indicTs = matrix(NA, nrow = max_designs, ncol = n)
	for (i_k in 1 : m){
		cat("  Kernel", i_k, "...")
		gd_res = resultsGreedySearch(initGreedyExperimentalDesignObject(
							X, max_designs = kernel_pre_num_designs, Kgram = Kgrams[[i_k]], objective = "kernel",
							num_cores = num_cores, seed = seed, start = TRUE, wait = TRUE, verbose = FALSE), 
						max_vectors = max_designs, form = "one_zero")
		cat(" done.\n")
		#and get the objective values for this kernel for all r designs
		# ending_indicTs from rJava is n x max_designs (inner arrays become columns); transpose to max_designs x n
		ending_indicTs_t = t(gd_res$ending_indicTs)
		objvalsi = compute_objective_vals_gpu(ending_indicTs_t, Kgrams[[i_k]])
		kernel_obj_values[, i_k] = objvalsi
		if (i_k == 1){ #start with the first kernel's best vectors
			initial_starting_indicTs = ending_indicTs_t
		}
	}
	
	#set the starting points
	for (r in 1 : max_designs){
		.jcall(java_obj, "V", "setStartingIndicT", as.integer(r - 1), as.integer(initial_starting_indicTs[r, ]))
	}
	
	#set the kernel weights
	.jcall(java_obj, "V", "setKernelWeights", kernel_weights)
	
	#set the maximum gain scaling
	if (maximum_gain_scaling){
		.jcall(java_obj, "V", "setMaximumGainScaling")
	}
	
	#feed in the gram matrices
	for (i_k in 1 : m){
		for (i in 1 : n){
			.jcall(java_obj, "V", "setSpecificKgramByRow", as.integer(i_k - 1), as.integer(i - 1), Kgrams[[i_k]][i, , drop = FALSE]) #java indexes from 0...n-1
		}
	}
	
	#feed in the initial objective values for each kernel
	for (r in 1 : max_designs){
		.jcall(java_obj, "V", "setInitialKernelObjValues", as.integer(r - 1), kernel_obj_values[r, ])
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
	greedy_multiple_kernel_experimental_design = list()
	greedy_multiple_kernel_experimental_design$max_designs = max_designs
	greedy_multiple_kernel_experimental_design$semigreedy = semigreedy
	greedy_multiple_kernel_experimental_design$start = start
	greedy_multiple_kernel_experimental_design$wait = wait
	greedy_multiple_kernel_experimental_design$diagnostics = diagnostics
	greedy_multiple_kernel_experimental_design$verbose = verbose
	greedy_multiple_kernel_experimental_design$X = X
	greedy_multiple_kernel_experimental_design$Kgrams = Kgrams
	greedy_multiple_kernel_experimental_design$n = n
	greedy_multiple_kernel_experimental_design$p = p
	greedy_multiple_kernel_experimental_design$m = m
	greedy_multiple_kernel_experimental_design$java_obj = java_obj
	class(greedy_multiple_kernel_experimental_design) = "greedy_multiple_kernel_experimental_design"
	
	#start if requested
	if (start){
		startSearch(greedy_multiple_kernel_experimental_design)
	}
	
	greedy_multiple_kernel_experimental_design
}

#' Returns the results of a greedy multiple kernel search
#' 
#' @param obj 				The \code{greedy_multiple_kernel_experimental_design} object where the search was run.
#' @param max_vectors		How many random allocation vectors you wish to return. The default is \code{NULL} indicating you want all of them.
#' @param form				Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's.
#' 
#' @author Adam Kapelner
#' @export
resultsMultipleKernelGreedySearch = function(obj, max_vectors = NULL, form = "one_zero"){
	assertClass(obj, "greedy_multiple_kernel_experimental_design")
	assertChoice(form, c("one_zero", "pos_one_min_one"))
	
	max_vectors_completed = greedySearchCurrentProgress(obj)
	if (is.null(max_vectors)){
		max_vectors = max_vectors_completed
	}
	assertInt(max_vectors, lower = 1, upper = max_vectors_completed)
	
	ending_indicTs = matrix(NA, nrow = max_vectors, ncol = obj$n)
	for (r in 1 : max_vectors){
		ending_indicTs[r, ] = .jcall(obj$java_obj, "[I", "getEndingIndicT", as.integer(r - 1))
	}
	
	if (form == "pos_one_min_one"){
		ending_indicTs = (ending_indicTs - 0.5) * 2
	}
	
	list(
		ending_indicTs = ending_indicTs
	)
}

#' Prints a summary of a \code{greedy_multiple_kernel_experimental_design} object
#'
#' @param object 	The \code{greedy_multiple_kernel_experimental_design} object to be summarized in the console
#' @param ...		Optional arguments for the \code{summary} function
#'
#' @author Adam Kapelner
#' @method summary greedy_multiple_kernel_experimental_design
#' @export
summary.greedy_multiple_kernel_experimental_design = function(object, ...){
	cat("Greedy Multiple Kernel Experimental Design Search\n")
	cat("  n =", object$n, "subjects\n")
	cat("  m =", object$m, "kernels\n")
	cat("  max_designs =", object$max_designs, "candidate designs\n")
	cat("  Current progress:", greedySearchCurrentProgress(object), "/", object$max_designs, "\n")
}

#' Prints a summary of a \code{greedy_multiple_kernel_experimental_design} object
#' 
#' @param x 		The \code{greedy_multiple_kernel_experimental_design} object to be printed in the console
#' @param ...		Optional arguments for the \code{print} function
#' 
#' @author Adam Kapelner
#' @method print greedy_multiple_kernel_experimental_design
#' @export
print.greedy_multiple_kernel_experimental_design = function(x, ...){
	summary(x, ...)
}

#' Plots a summary of a \code{greedy_multiple_kernel_experimental_design} object
#' 
#' @param x 		The \code{greedy_multiple_kernel_experimental_design} object to be plotted
#' @param ...		Optional arguments for the \code{plot} function
#' 
#' @author Adam Kapelner
#' @method plot greedy_multiple_kernel_experimental_design
#' @export
plot.greedy_multiple_kernel_experimental_design = function(x, ...){
	# Placeholder for plot implementation
	cat("Plotting for greedy_multiple_kernel_experimental_design objects not yet implemented.\n")
}
