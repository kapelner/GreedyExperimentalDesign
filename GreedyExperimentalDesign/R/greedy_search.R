#' Begin A Greedy Pair Switching Search
#' 
#' This method creates an object of type greedy_experimental_design and will immediately initiate
#' a search through allocation space for forced balance designs. For debugging, you can use set the \code{seed}
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
#' @param verbose			Should the algorithm emit progress output? Default is \code{TRUE}.
#' @param use_safe_inverse	Should a regularized inverse be used for the Mahalanobis objective?
#' 							Default is \code{FALSE}.
#' @return					An object of type \code{greedy_experimental_design_search} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 5,
#'   num_cores = 1,
#'   objective = "abs_sum_diff",
#'   start = TRUE,
#'   wait = TRUE,
#'   verbose = FALSE
#' )
#' ged
#' }
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
		seed = NULL,
		verbose = TRUE,
		use_safe_inverse = FALSE){
	
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
	assertLogical(verbose)
	assertLogical(use_safe_inverse)
	
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
			if (use_safe_inverse){
				SinvX = safe_cov_inverse(X)
			} else {
				SinvX = solve(stats::var(X))
			}
		}
	}
	

	
	#we are about to construct a GreedyExperimentalDesign java object. First, let R garbage collect
	#to clean up previous GreedyExperimentalDesign objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("GreedyExperimentalDesign.GreedyExperimentalDesign")
	set_verbose_if_available(java_obj, verbose)
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
	greedy_experimental_design_search$verbose = verbose
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
	progress = greedySearchCurrentProgress(x)
	res = resultsGreedySearch(x, max_vectors = 2)
	breaks = max(1, round(progress / 10))
	hist_data = data.frame(obj_val = res$obj_vals_orig_order)
	hist_plot = ggplot2::ggplot(hist_data, ggplot2::aes(x = obj_val)) +
		ggplot2::geom_histogram(bins = breaks, color = "black", fill = "grey70") +
		ggplot2::labs(x = "objective value", y = NULL, title = paste("After", progress, "searches")) +
		ggplot2::theme_classic()
#	hist(res$num_iters, br = progress / 10, xlab = "# of search iterations", ylab = NULL, main = "")
	
	#now do the plot of number of searches needed
	order_plot = build_obj_val_order_statistic_plot(x)
	draw_plot_pair(hist_plot, order_plot$plot)
	invisible(list(val_order_stats = order_plot$val_order_stats))
}

#' Plots the objective value by iteration
#' 
#' @param res 		Results from a greedy search object
#' @param runs 		A vector of run indices you would like to see plotted (default is to plot the first up to 9)
#' 
#' @author 			Adam Kapelner
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 5,
#'   num_cores = 1,
#'   diagnostics = TRUE,
#'   objective = "abs_sum_diff",
#'   start = TRUE,
#'   wait = TRUE,
#'   verbose = FALSE
#' )
#' res = resultsGreedySearch(ged, max_vectors = 2)
#' plot_obj_val_by_iter(res)
#' }
#' @export
plot_obj_val_by_iter = function(res, runs = NULL){
	if (is.null(res$obj_val_by_iters)){
		stop("You need to set diagnostics = TRUE on the search object.")
	}
	
	if (is.null(runs)){
		runs = 1 : min(length(res$obj_val_by_iters), 9)
	}
	num_to_plot = length(runs)
	
	plot_data = data.frame(run = character(0), iteration = integer(0), obj_val = numeric(0))
	for (run in runs){
		obj_vals = res$obj_val_by_iters[[run]]
		if (length(obj_vals) == 0){
			next
		}
		plot_data = rbind(plot_data, data.frame(
			run = paste("Run #", run),
			iteration = seq_along(obj_vals),
			obj_val = obj_vals
		))
	}
	plot_data$run = factor(plot_data$run, levels = paste("Run #", runs))
	if (nrow(plot_data) == 0){
		blank_plot = build_blank_plot("iteration", "objective value", title = NULL)
		print(blank_plot)
		return(invisible(NULL))
	}
	plot_data$line_ok = ave(plot_data$iteration, plot_data$run, FUN = function(x) length(x) > 1)
	iter_plot = ggplot2::ggplot(plot_data, ggplot2::aes(x = iteration, y = obj_val)) +
		ggplot2::geom_line(data = plot_data[plot_data$line_ok, , drop = FALSE]) +
		ggplot2::geom_point() +
		ggplot2::facet_wrap(~run, ncol = ceiling(sqrt(num_to_plot))) +
		ggplot2::labs(x = "iteration", y = "objective value") +
		ggplot2::theme_classic()
	print(iter_plot)
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
#' @examples
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 5,
#'   num_cores = 1,
#'   objective = "abs_sum_diff",
#'   start = TRUE,
#'   wait = TRUE,
#'   verbose = FALSE
#' )
#' plot_obj_val_order_statistic(ged, order_stat = 1, skip_every = 1)
#' }
#' @export
plot_obj_val_order_statistic = function(obj, order_stat = 1, skip_every = 5, type = "o", ...){
	order_plot = build_obj_val_order_statistic_plot(obj, order_stat = order_stat, skip_every = skip_every, type = type, ...)
	print(order_plot$plot)
	invisible(list(val_order_stats = order_plot$val_order_stats))	
}

build_blank_plot = function(xlab, ylab, title = NULL, xlim = NULL, ylim = NULL){
	blank_data = data.frame(x = 0, y = 0)
	blank_plot = ggplot2::ggplot(blank_data, ggplot2::aes(x = x, y = y)) +
		ggplot2::geom_blank() +
		ggplot2::labs(x = xlab, y = ylab, title = title) +
		ggplot2::theme_classic()
	if (!is.null(xlim) || !is.null(ylim)){
		blank_plot = blank_plot + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
	}
	blank_plot
}

draw_plot_pair = function(left_plot, right_plot){
	grid::grid.newpage()
	grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
	print(left_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
	print(right_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

build_obj_val_order_statistic_plot = function(obj, order_stat = 1, skip_every = 5, type = "o", ...){
	progress = greedySearchCurrentProgress(obj)
	ylab = paste("objective value (", order_stat, ")", sep = "")
	dots = list(...)
	title = NULL
	if (!is.null(dots$main)){
		title = dots$main
	}
	xlim = dots$xlim
	ylim = dots$ylim
	plot_color = NULL
	if (!is.null(dots$col)){
		plot_color = dots$col
	} else if (!is.null(dots$color)){
		plot_color = dots$color
	}
	line_size = dots$lwd
	point_shape = dots$pch
	point_size = NULL
	if (!is.null(dots$cex)){
		point_size = dots$cex * 1.5
	}
	if (progress <= 0 || progress < order_stat){
		blank_plot = build_blank_plot("Number of Searches", ylab, title = title, xlim = xlim, ylim = ylim)
		return(list(plot = blank_plot, val_order_stats = numeric(0)))
	}
	res = resultsGreedySearch(obj, max_vectors = 1)	#don't need many vectors; keep minimal for speed concerns
	vals = res$obj_vals_orig_order
	val_order_stats = array(NA, progress)
	skip_every = max(1, skip_every)
	indices = seq(order_stat, progress, by = skip_every)
	if (length(indices) == 0){
		indices = progress
	} else if (tail(indices, 1) != progress){
		indices = c(indices, progress)
	}
	for (d in indices){
		vals_d = vals[1 : d]
		vals_d = vals_d[is.finite(vals_d)]
		if (length(vals_d) == 0){
			next
		}
		if (order_stat == 1){
			val_order_stats[d] = min(vals_d)
		} else if (order_stat <= length(vals_d)){
			val_order_stats[d] = sort(vals_d)[order_stat]
		}
	}
	if (!any(is.finite(val_order_stats))){
		blank_plot = build_blank_plot("Number of Searches", ylab, title = title, xlim = xlim, ylim = ylim)
		return(list(plot = blank_plot, val_order_stats = val_order_stats))
	}
	plot_data = data.frame(searches = seq_len(progress), value = val_order_stats)
	plot_data = plot_data[is.finite(plot_data$value), , drop = FALSE]
	order_plot = ggplot2::ggplot(plot_data, ggplot2::aes(x = searches, y = value)) +
		ggplot2::labs(x = "Number of Searches", y = ylab, title = title) +
		ggplot2::theme_classic()
	type = tolower(type)
	draw_line = grepl("l", type) || grepl("o", type) || grepl("b", type) || grepl("c", type)
	draw_points = grepl("p", type) || grepl("o", type) || grepl("b", type)
	if (identical(type, "n")){
		draw_line = FALSE
		draw_points = FALSE
	}
	if (draw_line && nrow(plot_data) > 1){
		line_args = list(na.rm = TRUE)
		if (!is.null(plot_color)){
			line_args$color = plot_color
		}
		if (!is.null(line_size)){
			line_args$size = line_size
		}
		order_plot = order_plot + do.call(ggplot2::geom_line, line_args)
	}
	if (draw_points && nrow(plot_data) > 0){
		point_args = list(na.rm = TRUE)
		if (!is.null(plot_color)){
			point_args$color = plot_color
		}
		if (!is.null(point_size)){
			point_args$size = point_size
		}
		if (!is.null(point_shape)){
			point_args$shape = point_shape
		}
		order_plot = order_plot + do.call(ggplot2::geom_point, point_args)
	}
	if ((nrow(plot_data) == 0) || (!draw_line && !draw_points)){
		order_plot = order_plot + ggplot2::geom_blank()
	}
	if (!is.null(xlim) || !is.null(ylim)){
		order_plot = order_plot + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
	}
	list(plot = order_plot, val_order_stats = val_order_stats)	
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
#' \dontrun{
#' set.seed(1)
#' X = matrix(rnorm(20), nrow = 10)
#' ged = initGreedyExperimentalDesignObject(
#'   X,
#'   max_designs = 5,
#'   num_cores = 1,
#'   objective = "abs_sum_diff",
#'   start = TRUE,
#'   wait = TRUE,
#'   verbose = FALSE
#' )
#' res = resultsGreedySearch(ged, max_vectors = 2)
#' res$obj_vals
#' }
#' @export
resultsGreedySearch = function(obj, max_vectors = 9, form = "one_zero"){
	if (!is.null(max_vectors)){
		assertCount(max_vectors, positive = TRUE)
	}
	obj_vals = .jcall(obj$java_obj, "[D", "getObjectiveVals")
	num_iters = .jcall(obj$java_obj, "[I", "getNumIters")
	#these two are in order, so let's order the indicTs by the final objective values
	ordered_indices = order(obj_vals)
	num_completed = length(obj_vals)
	if (is.null(max_vectors)){
		last_index = num_completed
	} else {
		last_index = min(max_vectors, num_completed)
	}
	ending_indicTs = NULL
	starting_indicTs = NULL
	switches = NULL
	xbarj_diffs = NULL
	obj_val_by_iters = NULL
	pct_vec_same = NULL
	ordered_java_indices = integer(0)
	ordered_java_indices_j = NULL
	if (last_index > 0){
		ordered_java_indices = as.integer(ordered_indices[1 : last_index] - 1)
		ordered_java_indices_j = .jarray(ordered_java_indices)
		ending_indicTs = .jcall(obj$java_obj, "[[I", "getEndingIndicTs", ordered_java_indices_j, simplify = TRUE)
		if (form == "pos_one_min_one"){
			if (length(ending_indicTs) > 0 && min(ending_indicTs) >= 0 && max(ending_indicTs) <= 1){
				ending_indicTs = (ending_indicTs - 0.5) * 2
			}
		} else if (form == "one_zero" && length(ending_indicTs) > 0 && min(ending_indicTs) < 0){
			ending_indicTs = (ending_indicTs + 1) / 2
		}
	}
	if (obj$diagnostics && last_index > 0){
		starting_indicTs = .jcall(obj$java_obj, "[[I", "getStartingIndicTs", ordered_java_indices_j, simplify = TRUE)
		if (form == "pos_one_min_one"){
			if (length(starting_indicTs) > 0 && min(starting_indicTs) >= 0 && max(starting_indicTs) <= 1){
				starting_indicTs = (starting_indicTs - 0.5) * 2
			}
		} else if (form == "one_zero" && length(starting_indicTs) > 0 && min(starting_indicTs) < 0){
			starting_indicTs = (starting_indicTs + 1) / 2
		}
		switches = .jcall(obj$java_obj, "[[[I", "getSwitchedPairs", ordered_java_indices_j, simplify = TRUE)
		xbarj_diffs = .jcall(obj$java_obj, "[[[D", "getXbarjDiffs", ordered_java_indices_j, simplify = TRUE)
		obj_val_by_iters = .jcall(obj$java_obj, "[[D", "getObjValByIter", ordered_java_indices_j, simplify = TRUE)
		
		pct_vec_same = colSums(starting_indicTs == ending_indicTs) / length(starting_indicTs[, 1]) * 100
	}
	greedy_experimental_design_search_results = list(
		ordered_indices = ordered_indices,
		obj_vals = obj_vals[ordered_indices], 
		obj_vals_orig_order = obj_vals,
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
#' \dontrun{
#' optimize_asymmetric_treatment_assignment(n = 100)
#' optimize_asymmetric_treatment_assignment(n = 100, c_treatment = 2, c_control = 1)
#' optimize_asymmetric_treatment_assignment(c_total_max = 50, c_treatment = 2, c_control = 1)
#' }
#' @export
optimize_asymmetric_treatment_assignment = function(
		c_treatment = NULL,
		c_control = NULL,
		c_total_max = NULL,
		n = NULL
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
			nT = ceiling(n * c_treatment / (c_treatment + c_control))
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
