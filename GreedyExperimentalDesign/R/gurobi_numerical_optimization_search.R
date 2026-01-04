#' Begin Gurobi Optimized Search
#' 
#' This method creates an object of type optimal_experimental_design and will immediately initiate
#' a search through allocation space for forced balance designs. Make sure you setup Gurobi properly first. This means
#' applying for a license, downloading, installing, registering it on your computer using the \code{grbgetkey} command
#' with the license file in the default directory. Then, in R, install the package from the file in your gurobi directory. 
#' 
#' Currently, this method does not return multiple vectors. This will be improved in a later 
#' version. If you want this functionality now, use the hacked-up method \code{gurobi_multiple_designs}.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param objective			The objective function to use when searching design space. This is a string
#' 							with valid values "\code{mahal_dist}" (the default) or "\code{kernel}".
#' @param Kgram				If the \code{objective = kernel}, this argument is required to be an \code{n x n} matrix whose
#' 							entries are the evaluation of the kernel function between subject i and subject j. Default is \code{NULL}.
#' @param initial_time_limit_sec	The maximum amount of time the optimizer can run for in seconds. The default is \code{5 * 60}.
#' @param restart_time_limit_sec	The maximum amount of time each restart can run for in seconds. The default is \code{60}.
#' @param allow_restarts			Should the solver run additional restarts if too few unique solutions are returned?
#' 								Default is \code{TRUE}.
#' @param num_cores			Number of cores to use during search. Default is \code{2}.
#' @param w_0				The initial starting location (optional).
#' @param gurobi_params		A list of optional parameters to be passed to Gurobi (see their documentation online).
#' @param verbose			Should Gurobi log to console? Default is \code{TRUE}.
#' @param use_safe_inverse	Should a regularized inverse be used for the Mahalanobis objective?
#' 							Default is \code{FALSE}.
#' @param r					Number of solution vectors to request from the Gurobi pool.
#' @param pool_solutions 	Number of solutions to request from the Gurobi pool. Defaults to \code{10 * r}.
#' @param pool_gap 			Relative optimality gap for the pool. Default is \code{0.2}. Use \code{NULL} to skip.
#' @param pool_search_mode	Solution pool search mode. Default is \code{2} for diverse solutions.
#' @param mip_gap			Relative MIP gap target (stops when \code{|best-bound - best-incumbent| / |best-incumbent| <= mip_gap}).
#' 							Lower values force deeper search. Default is \code{1e-4}.
#' @param mip_gap_abs		Absolute MIP gap target (stops when \code{|best-bound - best-incumbent| <= mip_gap_abs}).
#' 							Lower values force deeper search. Default is \code{1e-10}.
#' @param mip_focus			Search focus: \code{0} (balance), \code{1} (find feasible solutions), \code{2} (prove optimality),
#' 							\code{3} (bound improvement). Default is \code{1}.
#' @param heuristics			Heuristics effort in \code{[0,1]} where higher values spend more time on heuristics.
#' 							Default is \code{0.2}.
#' @param cuts				Cut aggressiveness: \code{-1} (automatic), \code{0} (off), \code{1} (conservative),
#' 							\code{2} (aggressive), \code{3} (very aggressive). Default is \code{2}.
#' @param presolve			Presolve aggressiveness: \code{-1} (automatic), \code{0} (off), \code{1} (conservative),
#' 							\code{2} (aggressive). Default is \code{2}.
#' @return 					A list object which houses the results from Gurobi. Depending on the \code{gurobi_parms},
#' 							the data within will be different. The most relevant tags are \code{x} for the best found solution and \code{objval}
#' 							for the object
#' 
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' if ("gurobi" %in% loadedNamespaces()) {
#'   set.seed(1)
#'   X = matrix(rnorm(12), nrow = 6)
#'   gobj = initGurobiNumericalOptimizationExperimentalDesignObject(
#'     X,
#'     r = 2,
#'     num_cores = 1,
#'     initial_time_limit_sec = 5,
#'     verbose = FALSE
#'   )
#'   gobj$n
#' }
#' }
#' @export
initGurobiNumericalOptimizationExperimentalDesignObject = function(
		X = NULL, 
		objective = "mahal_dist",
		Kgram = NULL,
		num_cores = 2, 
		w_0 = NULL,
		initial_time_limit_sec = 5 * 60,
		restart_time_limit_sec = 60,
		allow_restarts = TRUE,
		verbose = TRUE,
		gurobi_params = list(),
		use_safe_inverse = FALSE,
		r,
		pool_solutions = NULL,
		pool_gap = 0.2,
		pool_search_mode = 2,
		mip_gap = 1e-4,
		mip_gap_abs = 1e-10,
		mip_focus = 1,
		heuristics = 0.2,
		cuts = 2,
		presolve = 2){
	
	gurobi_ns = get_gurobi_namespace()
	assertCount(num_cores, positive = TRUE)
	assertNumeric(initial_time_limit_sec, lower = 1)
	assertNumeric(restart_time_limit_sec, lower = 1)
	assertLogical(allow_restarts)
	assertLogical(verbose)
	assertClass(gurobi_params, "list")
	assertLogical(use_safe_inverse)
	assertCount(r, positive = TRUE)
	assertCount(pool_solutions, positive = TRUE, null.ok = TRUE)
	assertNumeric(pool_search_mode, lower = 0)
	assertNumeric(pool_gap, lower = 0, null.ok = TRUE)
	assertNumeric(mip_gap, lower = 0)
	assertNumeric(mip_gap_abs, lower = 0)
	assertNumeric(mip_focus, lower = 0, upper = 3)
	assertNumeric(heuristics, lower = 0, upper = 1)
	assertNumeric(cuts, lower = -1, upper = 3)
	assertNumeric(presolve, lower = -1, upper = 2)
	if (!is.null(w_0)){
		assertSetEqual(unique(w_0), c(1, -1))
		assertTRUE(sum(w_0) == 0)
	}
	
	gurobi_params$LogToConsole = as.numeric(verbose)
	gurobi_params$Threads = num_cores	
	gurobi_params$TimeLimit = initial_time_limit_sec
	if (is.null(gurobi_params$MIPGap)) {
		gurobi_params$MIPGap = mip_gap
	}
	if (is.null(gurobi_params$MIPGapAbs)) {
		gurobi_params$MIPGapAbs = mip_gap_abs
	}
	if (is.null(gurobi_params$MIPFocus)) {
		gurobi_params$MIPFocus = mip_focus
	}
	if (is.null(gurobi_params$Heuristics)) {
		gurobi_params$Heuristics = heuristics
	}
	if (is.null(gurobi_params$Cuts)) {
		gurobi_params$Cuts = cuts
	}
	if (is.null(gurobi_params$Presolve)) {
		gurobi_params$Presolve = presolve
	}
	if (is.null(pool_solutions)) {
		pool_solutions = 10 * r
	}
	if (!is.null(pool_solutions) && is.null(gurobi_params$PoolSolutions)) {
		gurobi_params$PoolSolutions = pool_solutions
	}
	if (!is.null(pool_solutions) && is.null(gurobi_params$PoolSearchMode)) {
		gurobi_params$PoolSearchMode = pool_search_mode
	}
	if (!is.null(pool_solutions) && !is.null(pool_gap) && is.null(gurobi_params$PoolGap)) {
		gurobi_params$PoolGap = pool_gap
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
	
	# Let Q = X^T SinvX X =or= K
	if (objective == "mahal_dist"){
		if (p < n){
			if (use_safe_inverse){
				SinvX = safe_cov_inverse(X)
			} else {
				SinvX = solve(stats::var(X))
			}
		}
		Q = X %*% SinvX %*% t(X)
	} else if (!is.null(Kgram)){
		Q = Kgram
	}
	one_vec = rep(1, n)
	
	#   minimize_w w^T Q w = (2b - 1)^T Q (2b - 1) = 4b^T Q b - 2 1^T Q b - 2 b^T Q 1 + 1^T Q 1 = b^T 4Q b - b^T 4Q 1 + 1^T Q 1 
	# = minimize_b b^T 4Q b - b^T 4Q 1 + 1^T Q 1 
	# minimize
	#        [b_1 b_2 ... b_n]^T Q [b_1 b_2 ... b_n] - 1^T Q [b_1 b_2 ... b_n]
	# subject to
	#        b_1 + b_2 + ... + b_n = n / 2,
	#        b_1, b_2, ..., b_n binary
	
	model 		= list()
	#quadratic component
	model$Q     = 8 * Q
	model$obj   = as.numeric(-4 * (Q %*% one_vec))
	model$ObjCon = as.numeric(t(as.matrix(one_vec)) %*% Q %*% as.matrix(one_vec)) #this will allow the objective value to be interpretable during verbose printouts
	#constraint component
	model$A     = t(as.matrix(one_vec))
	model$rhs   = n / 2
	model$sense = "="
	model$vtype = "B"
	if (!is.null(w_0)){
		model$start = w_0
	}

	#now return information as an object (just a list)
	gurobi_numerical_optimization_experimental_design_search = list()
	gurobi_numerical_optimization_experimental_design_search$X = X
	gurobi_numerical_optimization_experimental_design_search$SinvX = SinvX
	gurobi_numerical_optimization_experimental_design_search$n = n
	gurobi_numerical_optimization_experimental_design_search$p = p
	gurobi_numerical_optimization_experimental_design_search$Kgram = Kgram
	gurobi_numerical_optimization_experimental_design_search$objective = objective
	gurobi_numerical_optimization_experimental_design_search$w_0 = w_0
	gurobi_numerical_optimization_experimental_design_search$gurobi_params = gurobi_params
	gurobi_numerical_optimization_experimental_design_search$use_safe_inverse = use_safe_inverse
	gurobi_numerical_optimization_experimental_design_search$r = r
	gurobi_numerical_optimization_experimental_design_search$pool_solutions = pool_solutions
	gurobi_numerical_optimization_experimental_design_search$pool_gap = pool_gap
	gurobi_numerical_optimization_experimental_design_search$pool_search_mode = pool_search_mode
	gurobi_numerical_optimization_experimental_design_search$restart_time_limit_sec = restart_time_limit_sec
	gurobi_numerical_optimization_experimental_design_search$allow_restarts = allow_restarts
	gurobi_numerical_optimization_experimental_design_search$model = model
	class(gurobi_numerical_optimization_experimental_design_search) = "gurobi_numerical_optimization_experimental_design_search"
	#run the optimization and return the final object	
	return(c(gurobi_numerical_optimization_experimental_design_search, gurobi_call(gurobi_ns, model, gurobi_params)))

# 	#we are about to construct a GurobiNumericalOptimizeExperimentalDesign java object. First, let R garbage collect
# 	#to clean up previous objects that are no longer in use. This is important
# 	#because R's garbage collection system does not "see" the size of Java objects. Thus,
# 	#you are at risk of running out of memory without this invocation. 
# 	gc() #Delete at your own risk!	
	
# 	#now go ahead and create the Java object and set its information
# 	error_obj = NULL
# 	java_obj = .jnew("GurobiNumericalOptimizeExperimentalDesign.GurobiNumericalOptimizeExperimentalDesign")
# 	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))
# 	.jcall(java_obj, "V", "setN", as.integer(n))
# 	if (objective != "kernel"){
# 		p = ncol(X)
# 		.jcall(java_obj, "V", "setP", as.integer(p))
# 	}
# #	cat("time limit min: ", as.numeric(time_limit_min), "\n")
# 	.jcall(java_obj, "V", "setObjective", objective)
# 	.jcall(java_obj, "V", "setTimeLimitMin", as.numeric(time_limit_sec / 60))
# #	if (!is.null(node_limit)){
# #		if (node_limit <= 1){
# #			stop("Node limit must be one or more.")
# #		}
# #		.jcall(java_obj, "V", "setNodeLimit", as.numeric(round(node_limit)))
# #	}	
# #	if (!is.null(max_solutions)){
# #		if (max_solutions < 1){
# #			stop("Max solutions must be one or more.")
# #		}
# #		.jcall(java_obj, "V", "setMaxSolutions", as.integer(round(max_solutions)))
# #	}	
	
# 	if (!verbose){
# 		.jcall(java_obj, "V", "turnGurobiLogOff")
# 	}
# #	.jcall(java_obj, "V", "setLogFilename", log_file)
	
# 	#feed in the gram matrix if applicable
# 	if (!is.null(Kgram)){
# 		setGramMatrix(java_obj, Kgram)
# 	} else {
# 		#feed in the raw data
# 		for (i in 1 : n){	
# 			.jcall(java_obj, "V", "setDataRow", as.integer(i - 1), X[i, , drop = FALSE]) #java indexes from 0...n-1
# 		}
		
# 		#feed in the inverse var-cov matrix
# 		for (j in 1 : p){
# 			.jcall(java_obj, "V", "setInvVarCovRow", as.integer(j - 1), SinvX[j, , drop = FALSE]) #java indexes from 0...n-1
# 		}
# 	}
	
# 	#now return information as an object (just a list)
# 	gurobi_numerical_optimization_experimental_design_search = list()
# 	gurobi_numerical_optimization_experimental_design_search$X = X
# 	gurobi_numerical_optimization_experimental_design_search$SinvX = SinvX
# 	gurobi_numerical_optimization_experimental_design_search$n = n
# 	gurobi_numerical_optimization_experimental_design_search$p = p
# 	gurobi_numerical_optimization_experimental_design_search$Kgram = Kgram
# 	gurobi_numerical_optimization_experimental_design_search$objective = objective
# #	gurobi_numerical_optimization_experimental_design_search$max_solutions = max_solutions
# 	gurobi_numerical_optimization_experimental_design_search$time_limit_sec = time_limit_sec
# 	gurobi_numerical_optimization_experimental_design_search$java_obj = java_obj
# 	class(gurobi_numerical_optimization_experimental_design_search) = "gurobi_numerical_optimization_experimental_design_search"
	
# 	#now start search
# 	startSearch(gurobi_numerical_optimization_experimental_design_search)
	
# 	#return the final object
# 	gurobi_numerical_optimization_experimental_design_search
}




# #' Find multiple designs using Gurobi
# #' 
# #' This method searches through allocation space using Gurobi's optimization many times.
# #' It finds many different solutions by permuting the rows of the design matrix and 
# #' rerunning the optimization.
# #' 
# #' @param X 		The design matrix with $n$ rows (one for each subject) and $p$ columns 
# #' 					(one for each measurement on the subject). This is the design matrix you wish 
# #' 					to search for a more optimal design.
# #' @param r 				The number of vectors that should be returned.
# #' @param pool_solutions 	Number of solutions to request from the Gurobi pool. Defaults to \code{r}.
# #' @param pool_gap 		Relative optimality gap for the pool. Default is \code{0.1}. Use \code{NULL} to skip.
# #' @param pool_search_mode	Solution pool search mode. Default is \code{2} for diverse solutions.
# #' @param diverse			Should the returned solutions be greedily diversified by Hamming distance? Default is \code{TRUE}.
# #' @param ... 				Additional arguments to be passed to \code{initGurobiNumericalOptimizationExperimentalDesignObject}.
# #' @return			A matrix of allocation vectors of dimension \code{r x n}.	 
# #' 
# #' @author Kapelner
# #' @examples
# #' \dontrun{
# #' if ("gurobi" %in% loadedNamespaces()) {
# #'   set.seed(1)
# #'   X = matrix(rnorm(12), nrow = 6)
# #'   indicTs = gurobi_multiple_designs(
# #'     X,
# #'     r = 2,
# #'     num_cores = 1,
# #'     initial_time_limit_sec = 5,
# #'     verbose = FALSE
# #'   )
# #'   dim(indicTs)
# #' }
# #' }
# #' @export
# gurobi_multiple_designs = function(X, r, pool_solutions = NULL, pool_gap = 0.1, pool_search_mode = 2, diverse = TRUE, ...){
# 	n = nrow(X)
# 	assertCount(r, positive = TRUE)
# 	assertLogical(diverse)
# 	if (is.null(pool_solutions)){
# 		pool_solutions = r
# 	}
# 	assertCount(pool_solutions, positive = TRUE)
# 	assertNumeric(pool_search_mode, lower = 0)
# 	assertNumeric(pool_gap, lower = 0, null.ok = TRUE)

# 	select_diverse_rows = function(mat, k){
# 		if (nrow(mat) <= k){
# 			return(seq_len(nrow(mat)))
# 		}
# 		selected = 1L
# 		remaining = setdiff(seq_len(nrow(mat)), selected)
# 		while (length(selected) < k){
# 			min_dist = vapply(remaining, function(i){
# 				min(rowSums(abs(mat[i, , drop = FALSE] - mat[selected, , drop = FALSE])))
# 			}, numeric(1))
# 			pick = remaining[which.max(min_dist)]
# 			selected = c(selected, pick)
# 			remaining = setdiff(remaining, pick)
# 		}
# 		selected
# 	}

# 	dots = list(...)
# 	gnoed = do.call(
# 		initGurobiNumericalOptimizationExperimentalDesignObject,
# 		c(
# 			list(
# 				X,
# 				r = r,
# 				pool_solutions = pool_solutions,
# 				pool_gap = pool_gap,
# 				pool_search_mode = pool_search_mode
# 			),
# 			dots
# 		)
# 	)
# 	res = resultsGurobiNumericalOptimizeSearch(gnoed)
# 	indicTs = res$indicTs
# 	if (is.null(dim(indicTs))){
# 		indicTs = matrix(indicTs, nrow = 1)
# 	}

# 	if (nrow(indicTs) >= r){
# 		if (diverse && r > 1){
# 			idx = select_diverse_rows(indicTs, r)
# 			indicTs = indicTs[idx, , drop = FALSE]
# 		} else {
# 			indicTs = indicTs[seq_len(r), , drop = FALSE]
# 		}
# 	} else {
# 		warning("Only ", nrow(indicTs), " unique Gurobi pool solutions available; returning with replacement to reach ", r, ".")
# 		if (nrow(indicTs) == 0){
# 			stop("No solutions available from Gurobi.")
# 		}
# 		indicTs = indicTs[sample(seq_len(nrow(indicTs)), r, replace = TRUE), , drop = FALSE]
# 	}

# 	indicTs
# }

# #' Find multiple designs using Gurobi and returns the minimum
# #' 
# #' This method searches through allocation space using Gurobi's optimization many times.
# #' It finds many different solutions by permuting the rows of the design matrix and 
# #' rerunning the optimization.
# #' 
# #' @param X 		The design matrix with $n$ rows (one for each subject) and $p$ columns 
# #' 					(one for each measurement on the subject). This is the design matrix you wish 
# #' 					to search for a more optimal design.
# #' @param r 		The number of vectors that should be returned
# #' @param objective The objective function to calculate. Default is \code{"mahal_dist"} for Mahalanobis distance
# #' @param ... 		Additional arguments to be passed to \code{initGurobiNumericalOptimizationExperimentalDesignObject}.
# #' @return			A list with the minimum objective value and its vector	 
# #' 
# #' @author Kapelner
# #' @examples
# #' \dontrun{
# #' if ("gurobi" %in% loadedNamespaces()) {
# #'   set.seed(1)
# #'   X = matrix(rnorm(12), nrow = 6)
# #'   res = gurobi_min_of_multiple_designs(
# #'     X,
# #'     r = 2,
# #'     num_cores = 1,
# #'     initial_time_limit_sec = 5,
# #'     verbose = FALSE
# #'   )
# #'   res$obj
# #' }
# #' }
# #' @export
# gurobi_min_of_multiple_designs = function(X, r, objective = "mahal_dist", ...){
# 	indicTs = gurobi_multiple_designs(X, r, ...)
# 	obj_min = .Machine$double.xmax
# 	i_min = NA
# 	for (i in 1 : r){
# 		obj_i = compute_objective_val(X, indicTs[i, ], objective = "mahal_dist")
# 		if (obj_i < obj_min){
# 			obj_min = obj_i
# 			i_min = i
# 		}
# 	}
# 	list(indicT = indicTs[i_min, ], obj = obj_min)
# }

#' Query the Gurobi Results
#' 
#' Returns the results (thus far) of the Gurobi numerical optimization design search
#' 
#' @param obj 				The \code{gurobi_numerical_optimization_experimental_design_search} object that is currently running the search
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' if ("gurobi" %in% loadedNamespaces()) {
#'   set.seed(1)
#'   X = matrix(rnorm(12), nrow = 6)
#'   gobj = initGurobiNumericalOptimizationExperimentalDesignObject(
#'     X,
#'     r = 2,
#'     num_cores = 1,
#'     initial_time_limit_sec = 5,
#'     verbose = FALSE
#'   )
#'   res = resultsGurobiNumericalOptimizeSearch(gobj)
#'   res$obj_vals
#' }
#' }
#' @export
resultsGurobiNumericalOptimizeSearch = function(obj){
	indicTs = NULL
	pool_vecs = list()
	extract_pool_vec = function(sol){
		if (is.null(sol)){
			return(NULL)
		}
		if (is.numeric(sol)){
			return(sol)
		}
		if (is.list(sol)){
			for (field in c("poolx", "poolnx", "x", "xn")){
				val = sol[[field]]
				if (!is.null(val)){
					return(val)
				}
			}
		}
		NULL
	}
	make_balanced_start = function(n){
		w = integer(n)
		w[sample.int(n, n / 2)] = 1L
		w
	}
	coerce_pool_vecs = function(pool_val, n, solcount = NULL){
		if (is.null(pool_val)){
			return(list())
		}
		if (is.matrix(pool_val) || is.data.frame(pool_val)){
			pool_mat = as.matrix(pool_val)
			if (nrow(pool_mat) == n){
				return(lapply(seq_len(ncol(pool_mat)), function(j) pool_mat[, j]))
			}
			if (ncol(pool_mat) == n){
				return(lapply(seq_len(nrow(pool_mat)), function(i) pool_mat[i, ]))
			}
			return(list(as.numeric(pool_mat)))
		}
		if (is.list(pool_val)){
			flat = list()
			for (sol in pool_val){
				val = extract_pool_vec(sol)
				if (is.null(val)){
					next
				}
				if (is.matrix(val) || is.data.frame(val)){
					mat = as.matrix(val)
					if (nrow(mat) == n){
						flat = c(flat, lapply(seq_len(ncol(mat)), function(j) mat[, j]))
					} else if (ncol(mat) == n){
						flat = c(flat, lapply(seq_len(nrow(mat)), function(i) mat[i, ]))
					} else {
						flat = c(flat, list(as.numeric(mat)))
					}
				} else {
					flat = c(flat, list(val))
				}
			}
			return(flat)
		}
		if (is.numeric(pool_val)){
			if (!is.null(solcount) && length(pool_val) == n * solcount){
				mat = matrix(pool_val, nrow = n)
				return(lapply(seq_len(ncol(mat)), function(j) mat[, j]))
			}
			if (length(pool_val) %% n == 0 && length(pool_val) > n){
				mat = matrix(pool_val, nrow = n)
				return(lapply(seq_len(ncol(mat)), function(j) mat[, j]))
			}
			return(list(pool_val))
		}
		list()
	}
	solcount = obj$solcount
	if (is.null(solcount) && !is.null(obj$poolobjval)) {
		solcount = length(obj$poolobjval)
	}
	for (field in c("pool", "poolx", "poolnx", "poolX", "poolnX")){
		pool_vecs = coerce_pool_vecs(obj[[field]], obj$n, solcount)
		if (length(pool_vecs) > 0){
			break
		}
	}
	if (length(pool_vecs) == 0){
		for (field in c("xn", "x")){
			pool_vecs = coerce_pool_vecs(obj[[field]], obj$n, 1)
			if (length(pool_vecs) > 0){
				break
			}
		}
	}
	if (length(pool_vecs) > 0){
		pool_vecs = Filter(function(x) !is.null(x) && length(x) > 0, pool_vecs)
		pool_vecs = lapply(pool_vecs, function(x) as.integer(as.numeric(x) > 0.5))
		indicTs = ged_helper_normalize_rows(pool_vecs, n = obj$n)
	}

	if (is.null(indicTs)){
		stop("No Gurobi solution found on object.")
	}

	if (isTRUE(obj$allow_restarts) && !is.null(obj$r) && nrow(indicTs) < obj$r && !is.null(obj$model)) {
		gurobi_ns = tryCatch(get_gurobi_namespace(), error = function(e) NULL)
		if (!is.null(gurobi_ns)) {
			restarts = max(0L, obj$r - nrow(indicTs))
			base_seed = if (!is.null(obj$gurobi_params$Seed)) {
				as.integer(obj$gurobi_params$Seed)
			} else {
				sample.int(1000000L, 1L)
			}
			for (attempt in seq_len(restarts)) {
				model = obj$model
				model$start = make_balanced_start(obj$n)
				jitter_scale = 1e-6 * max(1, mean(abs(model$obj)))
				model$obj = as.numeric(model$obj) + stats::runif(length(model$obj), -jitter_scale, jitter_scale)
				params = obj$gurobi_params
				if (!is.null(obj$restart_time_limit_sec)) {
					params$TimeLimit = obj$restart_time_limit_sec
				}
				params$Seed = base_seed + attempt
				sol = gurobi_call(gurobi_ns, model, params)
				solcount = sol$solcount
				if (is.null(solcount) && !is.null(sol$poolobjval)) {
					solcount = length(sol$poolobjval)
				}
				new_vecs = list()
				for (field in c("pool", "poolx", "poolnx", "poolX", "poolnX")){
					new_vecs = coerce_pool_vecs(sol[[field]], obj$n, solcount)
					if (length(new_vecs) > 0){
						break
					}
				}
				if (length(new_vecs) == 0){
					for (field in c("xn", "x")){
						new_vecs = coerce_pool_vecs(sol[[field]], obj$n, 1)
						if (length(new_vecs) > 0){
							break
						}
					}
				}
				if (length(new_vecs) > 0){
					new_vecs = lapply(new_vecs, function(x) as.integer(as.numeric(x) > 0.5))
					new_mat = ged_helper_normalize_rows(new_vecs, n = obj$n)
					indicTs = rbind(indicTs, new_mat)
				}
				if (nrow(indicTs) >= obj$r){
					break
				}
			}
		}
	}

	indicTs = ged_helper_unique_rows(indicTs)
	if (nrow(indicTs) == 0){
		stop("No feasible Gurobi solutions found.")
	}
	message("Gurobi unique designs after restarts: ", nrow(indicTs))

	is_feasible = rowSums(indicTs) == obj$n / 2
	if (!all(is_feasible)) {
	  warning("Infeasible solutions violating n_T - n_C returned by Gurobi. Discarding.")
	  indicTs = indicTs[is_feasible, , drop = FALSE]
	}

	inv_cov_X = NULL
	if (identical(obj$objective, "mahal_dist") && !is.null(obj$SinvX)) {
		inv_cov_X = obj$SinvX
	}
	obj_vals = apply(indicTs, 1, FUN = function(wj){
		compute_objective_val(
			obj$X,
			wj,
			objective = obj$objective,
			inv_cov_X = inv_cov_X,
			use_safe_inverse = isTRUE(obj$use_safe_inverse)
		)
	})

	if (identical(obj$objective, "mahal_dist") && !is.null(obj$x)) {
		w_best = as.integer(as.numeric(obj$x) > 0.5)
		if (length(w_best) == obj$n) {
			objval_check = compute_objective_val(
				obj$X,
				w_best,
				objective = "mahal_dist",
				inv_cov_X = inv_cov_X,
				use_safe_inverse = isTRUE(obj$use_safe_inverse)
			)
			Q_check = if (!is.null(obj$Kgram)) {
				obj$Kgram
			} else {
				if (is.null(inv_cov_X)) {
					inv_cov_X = if (isTRUE(obj$use_safe_inverse)) {
						safe_cov_inverse(obj$X)
					} else {
						solve(stats::var(obj$X))
					}
				}
				obj$X %*% inv_cov_X %*% t(obj$X)
			}
			one_vec = rep(1, obj$n)
			qb = as.numeric(Q_check %*% w_best)
			q1 = as.numeric(Q_check %*% one_vec)
			obj_no_const = as.numeric(4 * crossprod(w_best, qb) - 4 * crossprod(w_best, q1))
			obj_with_const = obj_no_const + as.numeric(crossprod(one_vec, q1))
			obj_expected = (4 / (obj$n ^ 2)) * obj_with_const
			delta_expected = abs(objval_check - obj_expected)
			tol_expected = 1e-8 * max(1, abs(obj_expected))
			if (is.finite(delta_expected) && delta_expected > tol_expected) {
				warning(
					"Computed objective value differs from expected scaled Gurobi objective (",
					format(objval_check, digits = 6),
					" vs ",
					format(obj_expected, digits = 6),
					"). Check objective alignment and scaling."
				)
			}
		}
	}
	order_idx = order(obj_vals)
	indicTs = indicTs[order_idx, , drop = FALSE]
	list(indicT = indicTs[1, ], indicTs = indicTs, obj_vals = obj_vals[order_idx])	
}

ged_helper_unique_rows = function(mat){
	mat[!duplicated(as.data.frame(mat)), , drop = FALSE]
}

ged_helper_normalize_rows = function(mat, n = NULL){
	if (is.null(mat)){
		return(mat)
	}
	if (is.list(mat)){
		mat = do.call(rbind, lapply(mat, as.integer))
	} else if (!is.matrix(mat)){
		mat = matrix(mat, nrow = 1)
	}
	if (!is.null(n) && ncol(mat) != n && nrow(mat) == n){
		mat = t(mat)
	}
	if (is.null(dim(mat))){
		mat = matrix(mat, nrow = 1)
	}
	mat
}

get_gurobi_namespace = function() {
	if (!"gurobi" %in% loadedNamespaces()) {
		stop("Package \"gurobi\" is required. Call library(gurobi) before using this function.")
	}
	getNamespace("gurobi")
}

gurobi_call = function(gurobi_ns, model, params) {
	gurobi_fun = get("gurobi", envir = gurobi_ns)
	gurobi_fun(model, params)
}
