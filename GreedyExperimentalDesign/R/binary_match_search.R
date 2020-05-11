#' Begin a Search for Binary Matching Designs
#' 
#' This method creates an object of type binary_experimental_design and will find pairs. You can then
#' use the function \code{resultsBinaryMatchSearch} to create randomized allocation vectors. For one column
#' in X, we just sort to find the pairs trivially.
#' 
#' @param X					The design matrix with $n$ rows (one for each subject) and $p$ columns 
#' 							(one for each measurement on the subject). This is the design matrix you wish 
#' 							to search for a more optimal design.
#' @param compute_dist_matrix	The function that computes the distance matrix between every two observations in \code{X}, 
#' 								its only argument. The default is \code{NULL} signifying euclidean squared distance optimized in C++.
#' @param objective			The objective function to use when searching design space. This is a string
#' 							with valid values "\code{mahal_dist}" (the default), "\code{abs_sum_diff}".
#' @return					An object of type \code{binary_experimental_design} which can be further operated upon.
#' 
#' @author Adam Kapelner
#' @export
binaryMatchExperimentalDesignSearch = function(
		X,
		compute_dist_matrix = NULL,
		objective = "mahal_dist"){
	
	verify_objective_function(objective)
	
	n = nrow(X)
	p = ncol(X)
	if (n %% 2 != 0){
		stop("Design matrix must have even rows to have equal treatments and controls")
	}	
	if (objective == "mahal_dist" && p > n){
		stop("Objective cannot be mahalanobis distance if p > n.")
	}
	
	if (is.null(compute_dist_matrix) & p == 1){
		#we don't need to do anything except order them up
		indices_pairs = matrix(order(X[, 1]), ncol = 2, byrow = TRUE)
	} else {
		if (is.null(compute_dist_matrix)) {			
			euclidean_distance_sqd_cpp = CPP_FUNCTIONS[["euclidean_distance_sqd_cpp"]]
			if (is.null(euclidean_distance_sqd_cpp)){
				euclidean_distance_sqd_cpp = cppFunction('
					NumericMatrix compute_distance_matrix_cpp(NumericMatrix X) {
						int n = X.nrow();
						int p = X.ncol();
						NumericMatrix D(n, n);
						std::fill(D.begin(), D.end(), NA_REAL);
						for (int i_1 = 0; i_1 < (n - 1); i_1++){
							//Rcout << "computing for row #: " << (i_1 + 1) << "\\n";
							for (int i_2 = i_1 + 1; i_2 < n; i_2++){
								double sqd_diff = 0;
								for (int j = 0; j < p; j++){
									sqd_diff += pow(X(i_1, j) - X(i_2, j), 2); //by default the cmath library in std is loaded
								}
								D(i_1, i_2) = sqd_diff;
								D(i_2, i_1) = D(i_1, i_2);
							}
						}
						return D;
					}
				')
				CPP_FUNCTIONS[["euclidean_distance_sqd_cpp"]] = euclidean_distance_sqd_cpp #compile it once and cache
			}
			
			D = euclidean_distance_sqd_cpp(X)
		} else {
			D = compute_dist_matrix(X)
		}
		#ensure diagonal is infinity
		diag(D) = .Machine$double.xmax
		
		#get the matching solution using the heuristic
		matrix_sol = round(lp.assign(D)$solution) #just like Gurobi -- the 0/1's are sometimes not precisely 0/1
		
		indices_pairs = which(matrix_sol == 1, arr.ind = TRUE)
		for (i in 1 : n){
			indices_pairs[i, ] = sort(indices_pairs[i, ])
		}	
		indices_pairs = unique(indices_pairs)
	}

	#now return information as an object (just a list)
	binary_experimental_design = list()
	binary_experimental_design$X = X
	binary_experimental_design$n = n
	binary_experimental_design$p = p
	binary_experimental_design$objective = objective
	binary_experimental_design$compute_dist_matrix = compute_dist_matrix
	binary_experimental_design$D = D
	binary_experimental_design$indices_pairs = indices_pairs
	class(binary_experimental_design) = "binary_experimental_design"
	binary_experimental_design
}

CPP_FUNCTIONS = list()

#euclidean_distance_cpp = cppFunction('
#	NumericMatrix compute_distance_matrix_cpp(NumericMatrix X) {
#		int n = X.nrow();
#		int p = X.ncol();
#		NumericMatrix D(n, n);
#		std::fill(D.begin(), D.end(), NA_REAL);
#		for (int i_1 = 0; i_1 < (n - 1); i_1++){
#			//Rcout << "computing for row #: " << (i_1 + 1) << "\\n";
#			for (int i_2 = i_1 + 1; i_2 < n; i_2++){
#				double sqd_diff = 0;
#				for (int j = 0; j < p; j++){
#					sqd_diff += pow(X(i_1, j) - X(i_2, j), 2); //by default the cmath library in std is loaded
#				}
#				D(i_1, i_2) = sqrt(sqd_diff); //by default the cmath library in std is loaded
#				//D(i_2, i_1) = D(i_1, i_2);
#			}
#		}
#		return D;
#	}
#')

#' Returns allocation vectors that are binary matched
#' 
#' @param obj 				The \code{binary_experimental_design} object where the pairs are computed.
#' @param num_vectors		How many random allocation vectors you wish to return. The default is 1000.
#' @param compute_obj_vals	Should we compute all the objective values for each allocation? Default is \code{FALSE}.
#' @param form				Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' 
#' @author Adam Kapelner
#' @export
resultsBinaryMatchSearch = function(obj, num_vectors = 1000, compute_obj_vals = FALSE, form = "zero_one"){
	#now that we have the pairs, we can randomize for as many vectors as we wish
	indicTs = matrix(NA, nrow = num_vectors, ncol = obj$n)
	one_minus_one = matrix(
		sample(c(
			rep(1, num_vectors / 2 * n / 2), 
			rep(-1, num_vectors / 2 * n / 2)
		)), 
		nrow = num_vectors)
	
	minus_half_plus_half = c(-.5, .5)
	for (r in 1 : num_vectors){
		for (i in 1 : (n / 2)){
			indicTs[r, obj$indices_pairs[i, ]] = minus_half_plus_half * one_minus_one[r, i] + 0.5
		}
	}
	
	obj_vals = NULL
	if (compute_obj_vals){
		if (obj$objective == "mahal_dist"){
			SinvX = solve(var(obj$X))
			obj_vals = apply(indicTs, 1, FUN = function(w){compute_objective_val(obj$X, w, objective = "mahal_dist", SinvX)})
		} else {
			obj_vals = apply(indicTs, 1, FUN = function(w){compute_objective_val(obj$X, w, objective = obj$objective)})	
		}	
	}

	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	list(
		obj_vals_unordered = obj_vals,
		indicTs = indicTs
	)
}

#' Prints a summary of a \code{binary_experimental_design} object
#' 
#' @param x			The \code{binary_experimental_design} object to be summarized in the console
#' @param ...		Other parameters to pass to the default print function
#' 
#' @author 			Adam Kapelner
#' @method print binary_experimental_design
#' @export
print.binary_experimental_design = function(x, ...){
	cat("The pairs have been computed. Now use the resultsBinaryMatchSearch to make allocations.\n")
}

#' Prints a summary of a \code{binary_experimental_design} object
#' 
#' @param object		The \code{binary_experimental_design} object to be summarized in the console
#' @param ...			Other parameters to pass to the default summary function
#' 
#' @author 				Adam Kapelner
#' @method summary optimal_experimental_design_search
#' @export
summary.binary_experimental_design = function(object, ...){
	print(object, ...)
}
