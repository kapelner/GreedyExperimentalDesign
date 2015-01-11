
#' This method creates an object of type greedy_experimental_design and will immediately initiate
#' a search through $\indic{T}$ space.
#' 
#' @param X					The dataset you wish to search. We STRONGLY recommend you standardize your matrix 
#' 							using the method \code{\link{stdize_design_matrix}}
#' @param max_designs 		The maximum number of designs to be returned. Default is 10,000. Make this large 
#' 							so you can search however long you wish as the search can be stopped at any time by
#' 							using the \code{\link{stopGreedySearch}} method 
#' @param num_cores 		The number of CPU cores you wish to use during the search
#' @return					An object of type \code{greedy_experimental_design} which can be further operated upon
#' 
#' @author Adam Kapelner
#' @export
newGreedyExperimentalDesignObject = function(X, max_designs = 10000, num_cores = 1){
	#get dimensions immediately
	n = nrow(X)
	p = ncol(X)
	
	#now go ahead and create the Java object and set its information
	java_greedy_experimental_design = .jnew("GreedyExperimentalDesign.GreedyExperimentalDesign")
	.jcall(java_greedy_experimental_design, "V", "setMaxDesigns", max_designs)
	.jcall(java_greedy_experimental_design, "V", "setNumCores", num_cores)
	
	#now feed into Java some starting points for the search since it's easier to do in R
	for (d in 1 : max_designs){
		.jcall(java_greedy_experimental_design, "V", "addDesignStartPoint", create_random_dummy_vec(n))
	}
	
	#now begin the search
	.jcall(java_greedy_experimental_design, "V", "beginSearch")
		
	#now return information as the object
	greedy_experimental_design = list()
	greedy_experimental_design$java_greedy_experimental_design = java_greedy_experimental_design
	class(greedy_experimental_design) = "greedy_experimental_design"
	greedy_experimental_design
}



#' Creates a random binary vector which codes an experimental design
#' 
#' @param n			How many subjects does the experiment have 
#' @param p_w 		What is the proportion of treatments?
#' @return 			An array of length $n$ with $n \times p_w$ number of 1's and the rest 0's in a random order
#' 
#' @author Adam Kapelner
create_random_dummy_vec = function(n, p_w = 0.5){
	indic_T_dummy_permutation_vector = c(rep(1, n * p_w), rep(0, n * p_w)) #there are n_T 1's followed by n_C 0's dictated by p_w
	sample(indic_T_dummy_permutation_vector)
}


#' Standardizes a design matrix by subtracting each feature's mean and then dividing by each
#' feature's standard deviation
#' 
#' @param X			The design matrix to be standardized
#' @return 			The design matrix standardized of the same dimension as the original
#' 
#' @author Adam Kapelner
#' @export
stdize_design_matrix = function(X){
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}

#for debugging purposes only
generate_stdzied_design_matrix = function(n = 50, p = 3, covariate_dist = "iid_std_normal"){
	if (covariate_dist == "iid_std_uniform"){
		X = matrix(runif(n * p), nrow = n, ncol = p)	
	} else if (covariate_dist == "iid_std_normal"){
		X = matrix(rnorm(n * p), nrow = n, ncol = p)
	}
	#now standardize the matrix to make things easier later
	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
}

#compute_objectives = function(X, indic_T, inv_cov_X = NULL){
#	#saves computation sometimes to pass it in
#	if (is.null(inv_cov_X)){
#		inv_cov_X = solve(var(X))
#	}
#	
#	X_T = X[indic_T == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
#	X_C = X[indic_T == 0, , drop = FALSE] #coerce as matrix in order to matrix multiply later
#	X_T_bar = colMeans(X_T)
#	X_C_bar = colMeans(X_C)
#	X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar)
#	
#	abs_obj = sum(abs(X_T_bar_minus_X_C_bar))
#	#do the main calculation (eq 6) and return it as a scalar
#	
#	mahal_obj = as.numeric(t(X_T_bar_minus_X_C_bar) %*% inv_cov_X %*% X_T_bar_minus_X_C_bar)	
#	
#	list(abs_obj = abs_obj, mahal_obj = mahal_obj)
#}

#generate_stdzied_design_matrix = function(n, p, covariate_dist = "iid_std_normal"){
#	if (covariate_dist == "iid_std_uniform"){
#		X = matrix(runif(n * p), nrow = n, ncol = p)	
#	} else if (covariate_dist == "iid_std_normal"){
#		X = matrix(rnorm(n * p), nrow = n, ncol = p)
#	}
#	#now standardize the matrix to make things easier later
#	apply(X, 2, function(xj){(xj - mean(xj)) / sd(xj)})	
#}


alg_4 = function(X, n, MAX_ITER = 1000){
	cat("  4")
	#pick a random starting place in the space 
	indic_T_start = create_random_dummy_vec(n)	
	
	#create a matrix that will contain the history of the evolution of our indicator vector
	indic_T_history = matrix(NA, nrow = MAX_ITER, ncol = n)
	indic_T_history[1, ] = indic_T_start #we begin at the user-specified starting value
	#create a history of prop_D_M values, beginning with the initial
	abs_obj_history = array(NA, MAX_ITER)	
	abs_obj_history[1] = compute_objectives(X, indic_T_start)$abs_obj
	
#	cat("iter_num: 1  prop_D_M_min: ", prop_D_M_history[1], "indic_T_start: \n")
#	print(indic_T_start)
	
	#now we go through greedily until we reach the local minimum
	#we start at iteration num 2 because the first is set for us
	for (iter_num in 2 : MAX_ITER){	
		cat(".")
		#operate on the previous indic_T from last iteration
		indic_T = indic_T_history[iter_num - 1, ]		
		#initialize a new indic_T which will store the indic_T with the lower prop_D_M value
		#keep this initialized to NULL so we have a conveninent way of knowing if we didn't find one
		#that fit the bill
		indic_T_min = NULL
		#initialize the target obj to beat as the previous iteration's obj
		abs_obj_min = abs_obj_history[iter_num - 1]
		for (i_T in which(indic_T == 1)){
			for (i_C in which(indic_T == 0)){
				#make a copy of the indic_T vector 
				indic_T_proposal = indic_T
				#now make the single switch
				indic_T_proposal[i_T] = 0
				indic_T_proposal[i_C] = 1
				
				#R is too slow to call the function... so do bad programming!
				X_T = X[indic_T_proposal == 1, , drop = FALSE] #coerce as matrix in order to matrix multiply later
				X_C = X[indic_T_proposal == 0, , drop = FALSE] #coerce as matrix in order to matrix multiply later
				X_T_bar = colMeans(X_T)
				X_C_bar = colMeans(X_C)
				X_T_bar_minus_X_C_bar = as.matrix(X_T_bar - X_C_bar)
				#do the main calculation (eq 6) and return it as a scalar
				abs_obj = sum(abs(X_T_bar_minus_X_C_bar))
				
				#if this is lower than the current running winner...
				if (abs_obj < abs_obj_min){
					#then record the current running winning distance
					abs_obj_min = abs_obj
					#and record this proposal as the current running winner
					indic_T_min = indic_T_proposal
				}
			}
		}
		
		#if we didn't find any new indic_T - that means we couldn't beat the previous and thus we're done with the greedy search
		if (is.null(indic_T_min)){
			break
		}
		#record the new indic_T and the new minimum prop_D_M if we didn't break
		indic_T_history[iter_num, ] = indic_T_min
		abs_obj_history[iter_num] = abs_obj_min			
		
	}	
	cat("\n")
#	indic_T = indic_T_history[iter_num - 1, ]
#	l = compute_objectives(X, indic_T)
#	l[[indic_T]] = indic_T
#	l	
	compute_objectives(X, indic_T_history[iter_num - 1, ])
}

alg_4_many = function(X, n, w, MAX_ITER = 1000){	
	#indic_Ts = matrix(NA, nrow = n, ncol = w)
	abs_objs = array(NA, w)
	mahal_objs = array(NA, w)
	for (i in 1 : w){
		res = alg_4(X, n, MAX_ITER)
		#indic_Ts[, i] = res$indic_T
		abs_objs[i] = res$abs_obj
		mahal_objs[i] = res$mahal_obj
	}
	#now return the mins
	list(abs_obj = min(abs_objs), mahal_obj = min(mahal_objs))
}

run_all_algs_sims = function(ns = c(50, 100, 140, 200, 300, 500, 700, 1000), 
		ps = c(1, 2, 5, 10), 
		num_rerandomization_samples_rubin = 1,#1e5, 
		covariate_dist = "iid_std_normal",
		Nsim = 50){
	all_sims = matrix(NA, nrow = 0, ncol = 12)
	colnames(all_sims) = c(
			"n",
			"p",
			"covs_dist",
			"abs_obj_avg_alg_1",		
			"abs_obj_avg_alg_2",				
			"abs_obj_avg_alg_3",				
			"abs_obj_avg_alg_4",		
			"mahal_obj_avg_alg_1",
			"mahal_obj_avg_alg_2",
			"mahal_obj_avg_alg_3",
			"mahal_obj_avg_alg_4",
			"alg_3_n_T_over_n"
	)
	
	for (n in ns){
		for (p in ps){
			
			alg_1_res_abs_obj = array(NA, Nsim)
			alg_1_res_mahal_obj = array(NA, Nsim)
			alg_2_res_abs_obj = array(NA, Nsim)
			alg_2_res_mahal_obj = array(NA, Nsim)	
			alg_3_res_abs_obj = array(NA, Nsim)
			alg_3_n_ratio = array(NA, Nsim)
			alg_3_res_mahal_obj = array(NA, Nsim)
			alg_4_res_abs_obj = array(NA, Nsim)
			alg_4_res_mahal_obj = array(NA, Nsim)
			
			for (nsim in 1 : Nsim){
				#we want to average over the distribution of X too
				X = generate_stdzied_design_matrix(n, p)
				
				
				cat("nsim =", nsim, "n =", n, "p =", p, "\n")
				
#				alg_1_res = alg_1(X, n, num_rerandomization_samples_rubin)
#				alg_1_res_abs_obj[nsim] = alg_1_res$abs_obj
#				alg_1_res_mahal_obj[nsim] = alg_1_res$mahal_obj
				
				if (p == 1){
					alg_2C_res = alg_2_p_one(X, n)
				} else {
#					alg_2_res = alg_2_p_more_one(X, n)
				}
				
#				alg_2_res_abs_obj[nsim] = alg_2_res$abs_obj
#				alg_2_res_mahal_obj[nsim] = alg_2_res$mahal_obj	
				
				alg_3_res = alg_3(X, n)
				alg_3_res_abs_obj[nsim] = alg_3_res$abs_obj
				alg_3_res_mahal_obj[nsim] = alg_3_res$mahal_obj	
				alg_3_n_ratio[nsim] = alg_3_res$n_T_over_n
				
#				alg_4_res = alg_4(X, n)
#				alg_4_res_abs_obj[nsim] = alg_4_res$abs_obj
#				alg_4_res_mahal_obj[nsim] = alg_4_res$mahal_obj
			}
			
			####tally all results in final matrix
			all_sims = rbind(all_sims, c(
							n, 
							p, 
							covariate_dist, 
							mean(alg_1_res_abs_obj),				
							mean(alg_2_res_abs_obj),				
							mean(alg_3_res_abs_obj),				
							mean(alg_4_res_abs_obj),
							mean(alg_1_res_mahal_obj),
							mean(alg_2_res_mahal_obj),
							mean(alg_3_res_mahal_obj),
							mean(alg_4_res_mahal_obj),
							mean(alg_3_n_ratio)
					))
			#and save to HD iteratively so there's not too much loss if it crashes midway
			write.csv(all_sims, paste("all_sims_", n, ".csv", sep = ""), row.names = FALSE)			
		}
	}
}

run_rerand_alg_vs_greedy_sims = function(ns = c(40, 100, 140, 200, 300, 500, 700, 1000), 
		ps = c(1, 2, 5, 10), 
		covariate_dist = "iid_std_normal", 
		rs = c(10, 100, 316, 1e3, 1e5, 1e6),
		ws = c(2, 5, 10, 50, 100, 200),
		Nsim = 50){
	
	all_sims = matrix(NA, nrow = 0, ncol = 9)
	colnames(all_sims) = c(
			"n",
			"p",
			"covs_dist",
			"r",
			"abs_obj_avg_alg_1",
			"mahal_obj_avg_alg_1",
			"w",						
			"abs_obj_avg_alg_4",			
			"mahal_obj_avg_alg_4"
	)
	
	if (length(rs) != length(ws)){
		stop("there should be equal rs as ws")
	}
	
	for (n in ns){
		for (p in ps){
			for (j in 1 : length(rs)){
				r = rs[j]
				w = ws[j]
				
				alg_1_res_abs_obj = array(NA, Nsim)
				alg_1_res_mahal_obj = array(NA, Nsim)
				alg_4_res_abs_obj = array(NA, Nsim)
				alg_4_res_mahal_obj = array(NA, Nsim)
				
				for (nsim in 1 : Nsim){
					#we want to average over the distribution of X too
					X = generate_stdzied_design_matrix(n, p)					
					
					cat("nsim =", nsim, "n =", n, "p =", p, "r =", r, "w =", w, "\n")
					
					alg_1_res = alg_1(X, n, r)
					alg_1_res_abs_obj[nsim] = alg_1_res$abs_obj
					alg_1_res_mahal_obj[nsim] = alg_1_res$mahal_obj				
					
					alg_4_res = alg_4_many(X, n, w)
					alg_4_res_abs_obj[nsim] = alg_4_res$abs_obj
					alg_4_res_mahal_obj[nsim] = alg_4_res$mahal_obj
				}
				
				####tally all results in final matrix
				all_sims = rbind(all_sims, c(
								n, 
								p, 
								covariate_dist, 
								r,					
								mean(alg_1_res_abs_obj),
								mean(alg_1_res_mahal_obj),	 
								w,			
								mean(alg_4_res_abs_obj),
								mean(alg_4_res_mahal_obj)
						))
				#and save to HD iteratively so there's not too much loss if it crashes midway
				write.csv(all_sims, paste("greedy_vs_rerand_sims_", n, ".csv", sep = ""), row.names = FALSE)	
			}
		}
	}
}