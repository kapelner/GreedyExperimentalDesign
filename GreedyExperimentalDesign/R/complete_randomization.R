#' Implements forced balanced randomization
#' 
#' For debugging, you can use \code{set.seed}
#' to be assured of deterministic output.
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization_with_forced_balanced = function(n, r, form = "one_zero"){
	indicTs = matrix(NA, nrow = r, ncol = n)
	zero_one_vec = c(rep(0, n / 2), rep(1, n / 2))
	for (nsim in 1 : r){
		indicTs[nsim, ] = shuffle_cpp(zero_one_vec)
	}
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}


#' Implements complete randomization (without forced balance)
#' 
#' For debugging, you can use \code{set.seed}
#' to be assured of deterministic output.
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
complete_randomization = function(n, r, form = "one_zero"){
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	for (nsim in 1 : r){
		indicTs[nsim, ] = rbinom(n, 1, 0.5)
	}
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}

#' Implements unequally allocated complete randomization
#' 
#' For debugging, you can use \code{set.seed}
#' to be assured of deterministic output.
#' 
#' @param n 		number of observations
#' @param prop_T    the proportion of treatments needed
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
imbalanced_complete_randomization = function(n, prop_T, r, form = "one_zero"){
	assert_count(n, positive = TRUE)
	assert_numeric(prop_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
	n_T = n * prop_T
	assert_count(n_T, positive = TRUE)
	assert_count(r, positive = TRUE)
	assert_choice(form, c("one_zero", "pos_one_min_one"))
	
	indicTs = matrix(NA, nrow = r, ncol = n)
	
	blank = c(rep(1, n_T), rep(0, n - n_T))
	for (nsim in 1 : r){
		indicTs[nsim, ] = sample(blank)
	}
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}

#' Implements unequally allocated block designs
#' 
#' For debugging, you can use \code{set.seed}
#' to be assured of deterministic output. The following quantities
#' in this design must be integer valued or an error will be thrown: 
#'   n_B := n / B and n_B * prop_T
#' 
#' @param n 		number of observations
#' @param prop_T    the proportion of treatments needed
#' @param B 		the number of blocks
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @export
imbalanced_block_designs = function(n, prop_T, B, r, form = "one_zero"){
	assertCount(n, positive = TRUE)
	assert_numeric(prop_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
	assertCount(B, positive = TRUE) 
	assertCount(r, positive = TRUE)
	
	n_B = n / B
	n_B_T = n_B * prop_T
	assertCount(n_B, positive = TRUE)
	assertCount(n_B_T, positive = TRUE)
	dummy_block = c(rep(1, n_B_T), rep(0, n_B - n_B_T))
	
	Ws = list()
	for (b in 1 : B){
		Ws[[b]] = matrix(NA, nrow = r, ncol = n_B)
		for (nr in 1 : r){
			Ws[[b]][nr, ] = sample(dummy_block)
		}
	}
	indicTs = list.cbind(Ws)
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}