#' Implements forced balanced randomization
#' 
#' For debugging, you can use \code{set.seed}
#' to be assured of deterministic output.
#' 
#' @param n 		number of observations
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @param seed		An integer which is the seed to be set within C++. Default is \code{NULL} which means the seed is set from the system clock.
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' complete_randomization_with_forced_balanced(n = 6, r = 2, seed = 1)
#' }
#' @export
complete_randomization_with_forced_balanced = function(n, r, form = "one_zero", seed = NULL){
	assert_count(n, positive = TRUE)
	assert_count(r, positive = TRUE)
	assert_choice(form, c("one_zero", "pos_one_min_one"))
	assert_count(seed, positive = TRUE, null.ok = TRUE)
	seed = ifelse(is.null(seed), NA_integer_, as.integer(seed))
	assert_integer(seed)
	
	indicTs = complete_randomization_forced_balanced_cpp(n, r, seed)
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
#' @examples
#' \dontrun{
#' complete_randomization(n = 6, r = 2)
#' }
#' @export
complete_randomization = function(n, r, form = "one_zero"){
	assert_count(n, positive = TRUE)
	assert_count(r, positive = TRUE)
	assert_choice(form, c("one_zero", "pos_one_min_one"))
	
	indicTs = complete_randomization_cpp(n, r)
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
#' @param seed      An integer which is the seed to be set within C++. Default is \code{NULL} which means the seed is set from the system clock.
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' imbalanced_complete_randomization(n = 10, prop_T = 0.3, r = 2, seed = 1)
#' }
#' @export
imbalanced_complete_randomization = function(n, prop_T, r, form = "one_zero", seed = NULL){
	assert_count(n, positive = TRUE)
	assert_numeric(prop_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
	n_T = n * prop_T
	assert_count(n_T, positive = TRUE)
	assert_count(r, positive = TRUE)
	assert_choice(form, c("one_zero", "pos_one_min_one"))
	seed = ifelse(is.null(seed), NA_integer_, as.integer(seed))
	assert_integer(seed)
	
	indicTs = complete_randomization_imbalanced_cpp(n, n_T, r, seed)
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
#' @param prop_T    the proportion of treatments allocated
#' @param B 		the number of blocks
#' @param r 		number of randomized designs you would like
#' @param form		Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's. 
#' @param seed      An integer which is the seed to be set within C++. Default is \code{NULL} which means the seed is set from the system clock.
#' @return 			a matrix where each column is one of the \code{r} designs
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' imbalanced_block_designs(n = 12, prop_T = 0.5, B = 3, r = 2, seed = 1)
#' }
#' @export
imbalanced_block_designs = function(n, prop_T, B, r, form = "one_zero", seed = NULL){
	assert_count(n, positive = TRUE)
	assert_numeric(prop_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
	assert_count(B, positive = TRUE) 
	assert_count(r, positive = TRUE)
	seed = ifelse(is.null(seed), NA_integer_, as.integer(seed))
	assert_integer(seed)
	
	n_B = n / B
	n_B_T = n_B * prop_T
	assert_count(n_B, positive = TRUE)
	assert_count(n_B_T, positive = TRUE)
	indicTs = imbalanced_block_designs_cpp(n_B, n_B_T, B, r, seed)
	if (form == "pos_one_min_one"){
		indicTs = (indicTs - 0.5) * 2
	}
	indicTs
}

# internal cache for block design var-cov results
.gen_var_cov_cache = new.env(parent = emptyenv())
.gen_var_cov_cache$max = 8L

#' Computes varcov matrix for block designs
#' 
#' The varcov matrix for block designs consists of a block-
#' diagonal matrix with B blocks (the number of blocks in
#' the design) with off-diagonal entries = -1 / (n/B - 1)
#' where n is the number of subjected in the study.
#'
#' @param n 		number of observations
#' @param prop_T    the proportion of treatments allocated
#' @param B 		the number of blocks
#' @param use_cache	Cache results for repeated calls with identical inputs. Default is \code{TRUE}.
#' @return 			varcov matrix for the specific block design
#' 
#' @author Adam Kapelner
#' @examples
#' \dontrun{
#' gen_var_cov_matrix_block_designs(n = 12, prop_T = 0.5, B = 3)
#' }
#' @export
gen_var_cov_matrix_block_designs = function(n, prop_T, B, use_cache = TRUE){
  assertCount(n, positive = TRUE)
  assert_numeric(prop_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
  assertCount(B, positive = TRUE)
  assertLogical(use_cache)
  key = NULL
  if (use_cache){
    key = paste("cache", n, sprintf("%.17g", prop_T), B, sep = "|")
    cached = .gen_var_cov_cache[[key]]
    if (!is.null(cached)){
      return(cached)
    }
  }
  SigmaW = gen_var_cov_matrix_block_designs_cpp(n, prop_T, B)
  if (use_cache){
    .gen_var_cov_cache[[key]] = SigmaW
    order = .gen_var_cov_cache$order
    if (is.null(order)){
      order = character()
    }
    order = c(order[order != key], key)
    if (length(order) > .gen_var_cov_cache$max){
      drop_key = order[1]
      rm(list = drop_key, envir = .gen_var_cov_cache)
      order = order[-1]
    }
    .gen_var_cov_cache$order = order
  }
  SigmaW
}
