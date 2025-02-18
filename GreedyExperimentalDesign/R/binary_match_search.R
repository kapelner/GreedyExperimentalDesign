#' Compute Binary Matching Strcuture
#' 
#' This method creates an object of type binary_match_structure and will compute pairs. You can then
#' use the functions \code{initBinaryMatchExperimentalDesignSearch} and \code{resultsBinaryMatchSearch} 
#' to create randomized allocation vectors. For one column in X, we just sort to find the pairs trivially.
#' 
#' @param X						The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 								(one for each measurement on the subject). This is the design matrix you wish 
#' 								to search for a more optimal design.
#' @param compute_dist_matrix	The function that computes the distance matrix between every two observations in \code{X}, 
#' 								its only argument. The default is \code{NULL} signifying euclidean squared distance optimized in C++.
#' @param mahal_match			Match using Mahalanobis distance. Default is \code{FALSE}.
#' @param D						A distance matrix precomputed. The default is \code{NULL} indicating the distance matrix should be computed.
#' @return						An object of type \code{binary_experimental_design} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
computeBinaryMatchStructure = function(X, mahal_match = FALSE, compute_dist_matrix = NULL, D = NULL){
	assertClass(X, "matrix")
	assertClass(compute_dist_matrix, "function", null.ok = TRUE)
	n = nrow(X)
	assertTRUE(n > 1)
	p = ncol(X)
	if (!is.null(D)){
		assertClass(D, "matrix")
		assertTRUE(nrow(D) == n)
		assertTRUE(ncol(D) == n)
		for (i in 1 : (n - 1)){ #ensure symmetric
			for (j in (i + 1) : n){
				assertTRUE(D[i, j] == D[j, i])
			}
		}
	}
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}
	
	if (is.null(compute_dist_matrix) & is.null(D) & p == 1){
		#we don't need to do anything except order them up
		indicies_pairs = matrix(order(X[, 1]), ncol = 2, byrow = TRUE)
	} else {
		if (is.null(compute_dist_matrix) & is.null(D)) {	
			if (mahal_match){
				#C++-optimize one day please!
				S_X_inv = solve(var(X))
				D = matrix(NA, nrow = n, ncol = n)
				for (i in 1 : (n - 1)){ #ensure symmetric
					for (j in (i + 1) : n){
						xdiff = X[i, ] - X[j, ]
						D[i, j] = D[j, i] = xdiff %*% (S_X_inv %*% xdiff)
					}
				}
			} else {
				#default is C++-optimized sqd euclidean distance function		
				D = compute_distance_matrix_cpp(X)				
			}
		} else if (is.null(D)){
			D = compute_dist_matrix(X)
		}
		#ensure diagonal is infinity in the R language
		diag(D) = .Machine$double.xmax
		
		#get the matching solution using the heuristic
		indicies_pairs = as.matrix(nbpMatching::nonbimatch(nbpMatching::distancematrix(D))$matches[, c("Group1.Row", "Group2.Row")])
		
		for (i in 1 : n){
			indicies_pairs[i, ] = sort(indicies_pairs[i, ])
		}	
		indicies_pairs = unique(indicies_pairs)
	}

	#now return information as an object (just a list)
	binary_match_structure = list()
	binary_match_structure$X = X
	binary_match_structure$n = n
	binary_match_structure$p = p
	binary_match_structure$compute_dist_matrix = compute_dist_matrix
	binary_match_structure$D = D
	binary_match_structure$indicies_pairs = indicies_pairs
	class(binary_match_structure) = "binary_match_structure"
	binary_match_structure
}

#' Begin a Binary Match Search
#' 
#' This method creates an object of type pairwise_matching_experimental_design_search and will immediately initiate
#' a search through $1_{T}$ space for pairwise match designs based on the structure computed in the function \code{computeBinaryMatchStructure}. 
#' For debugging, you can use set the \code{seed} parameter and \code{num_cores = 1} to be assured of deterministic output.
#' 
#' @param binary_match_structure 	The \code{binary_experimental_design} object where the pairs are computed.
#' @param max_designs				How many random allocation vectors you wish to return. The default is 1000.
#' @param wait						Should the \code{R} terminal hang until all \code{max_designs} vectors are found? The 
#' 									default is \code{FALSE}.
#' @param start						Should we start searching immediately (default is \code{TRUE}).
#' @param num_cores					The number of CPU cores you wish to use during the search. The default is \code{1}.
#' @param seed						The set to set for deterministic output. This should only be set if \code{num_cores = 1} otherwise
#' 									the output will not be deterministic. Default is \code{NULL} for no seed set. 
#' @param prop_flips				Proportion of flips. Default is all. Lower for more correlated assignments (useful for research only).
#' 
#' @author Adam Kapelner
#' @export
initBinaryMatchExperimentalDesignSearch = function(binary_match_structure, 
		max_designs = 1000, 
		wait = FALSE, 
		start = TRUE,
		num_cores = 1,
		seed = NULL,
		prop_flips = 1){
	assertClass(binary_match_structure, "binary_match_structure")
	assertCount(max_designs, positive = TRUE)
	assertNumeric(prop_flips, lower = 0, upper = 1)
	
	if (prop_flips < 1){
		warning("prop_flips feature is not implemented yet")
	}
		
	n = binary_match_structure$n
	if (2^(n / 2) < max_designs){
		stop(paste("The total number of unique pairwise matching vectors is only", 2^(n / 2), "which is less than the", max_designs, "you requested."))
	}
	
	if (max_designs > 2^(n / 2) / 2){
		warning(paste("The total number of unique pairwise matching vectors is only", 2^(n / 2), "which is less than double the", max_designs, "you requested.\nThis can be very slow as the space of designs gets crowded quickly and it is difficult to find new ones."))
	}
		
	#we are about to construct a PairwiseMatchingExperimentalDesign java object. First, let R garbage collect
	#to clean up previous java design objects that are no longer in use. This is important
	#because R's garbage collection system does not "see" the size of Java objects. Thus,
	#you are at risk of running out of memory without this invocation. 
	gc() #Delete at your own risk!	
	
	#now go ahead and create the Java object and set its information
	java_obj = .jnew("PairwiseMatchingExperimentalDesign.PairwiseMatchingExperimentalDesign")
	.jcall(java_obj, "V", "setMaxDesigns", as.integer(max_designs))
	.jcall(java_obj, "V", "setNumCores", as.integer(num_cores))	
	if (!is.null(seed)){
		.jcall(java_obj, "V", "setSeed", as.integer(seed))
		if (num_cores != 1){
			warning("Setting the seed with multiple cores does not guarantee deterministic output.")
		}		
	}
	.jcall(java_obj, "V", "setN", as.integer(n))
	if (wait){
		.jcall(java_obj, "V", "setWait")
	}
	
	#now upload the structure to the java object
	for (i in 1 : (n / 2)){
		.jcall(java_obj, "V", "setMatchPairIndicies", 
			as.integer(i - 1), 
			as.integer(binary_match_structure$indicies_pairs[i, ] - 1))		
	}
	
	#now return information as an object (just a list)
	pairwise_matching_experimental_design_search = list()
	pairwise_matching_experimental_design_search$binary_match_structure = binary_match_structure
	pairwise_matching_experimental_design_search$max_designs = max_designs
	pairwise_matching_experimental_design_search$start = start
	pairwise_matching_experimental_design_search$wait = wait
	pairwise_matching_experimental_design_search$num_cores = num_cores
	pairwise_matching_experimental_design_search$java_obj = java_obj
	class(pairwise_matching_experimental_design_search) = "pairwise_matching_experimental_design_search"
	#if the user wants to run it immediately...
	if (start){
		startSearch(pairwise_matching_experimental_design_search)
	}
	#return the final object
	pairwise_matching_experimental_design_search
}

#' Binary Pair Match Search
#' 
#' Returns the results (thus far) of the binary pair match design search
#' 
#' @param obj 					The \code{pairwise_matching_experimental_design_search} object that is currently running the search
#' @param form					Which form should the assignments be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @export
resultsBinaryMatchSearch = function(obj, form = "one_zero"){
	assertClass(obj, "pairwise_matching_experimental_design_search")
	assertChoice(form, c("zero_one", "pos_one_min_one"))
	
	ending_indicTs = .jcall(obj$java_obj, "[[I", "getEndingIndicTs", simplify = TRUE)	
	if (form == "pos_one_min_one"){
		ending_indicTs = (ending_indicTs - 0.5) * 2
	}	
	ending_indicTs
}

#' Prints a summary of a \code{pairwise_matching_experimental_design_search} object
#' 
#' @param x			The \code{pairwise_matching_experimental_design_search} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print 	pairwise_matching_experimental_design_search
#' @export
print.pairwise_matching_experimental_design_search = function(x, ...){
	progress = .jcall(x$java_obj, "I", "progress")
	time_elapsed = searchTimeElapsed(x)
	if (progress == 0){
		cat("No progress on the PairwiseMatchingExperimentalDesign. Did you run \"startSearch?\"\n")
	} else if (progress == x$max_designs){
		cat("The search completed in", time_elapsed, "seconds.", progress, "vectors have been found.\n")
	} else {
		cat("The search has found ", progress, " vectors thus far (", round(progress / x$max_designs * 100), "%) in ", time_elapsed," seconds.\n", sep = "")
	}
}

#' Prints a summary of a \code{pairwise_matching_experimental_design_search} object
#' 
#' @param object		The \code{pairwise_matching_experimental_design_search} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary pairwise_matching_experimental_design_search
#' @export
summary.pairwise_matching_experimental_design_search = function(object, ...){
	print(object, ...)
}

#' Prints a summary of a \code{binary_match_structure} object
#' 
#' @param x			The \code{binary_match_structure} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print binary_match_structure
#' @export
print.binary_match_structure = function(x, ...){
	cat("The pairs have been computed. Now use the initBinaryMatchExperimentalDesignSearch function to create allocations.\n")
}

#' Prints a summary of a \code{binary_match_structure} object
#' 
#' @param object		The \code{binary_match_structure} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary binary_match_structure
#' @export
summary.binary_match_structure = function(object, ...){
	print(object, ...)
}
