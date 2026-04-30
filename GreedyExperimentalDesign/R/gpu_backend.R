#' Check for optional GPU support
#'
#' @return A logical scalar. Attribute \code{wgpu_compiled} is \code{TRUE} when
#'   the package was compiled with WebGPU backend support.
#' @export
ged_gpu_available = function(){
	ged_gpu_available_cpp()
}

#' List optional GPU devices
#'
#' @return A data frame with columns \code{id}, \code{name}, \code{backend},
#'   \code{type}, and \code{usable}. It is empty when GPU support was not
#'   compiled into the package.
#' @export
ged_gpu_devices = function(){
	ged_gpu_devices_cpp()
}

#' Compute a squared Euclidean distance matrix through the optional native backend
#'
#' @param X A numeric matrix.
#' @param backend One of \code{"auto"}, \code{"cpu"}, or \code{"gpu"}.
#' @param device Zero-based device id. -1 for native C++.
#' @return An \code{n x n} squared Euclidean distance matrix.
#' @export
compute_distance_matrix_gpu = function(X, backend = "auto", device = 0){
	assertMatrix(X, mode = "numeric", any.missing = FALSE)
	assertChoice(backend, c("auto", "cpu", "gpu"))
	assertInt(device, lower = -1)
	if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE)) || backend == "cpu") {
		backend = "cpu"
	}
	compute_distance_matrix_gpu_cpp(X, backend, as.integer(device))
}

#' Compute a kernel matrix through the optional native backend
#'
#' @param X A numeric matrix.
#' @param kernel One of \code{"gaussian"}, \code{"laplacian"},
#'   \code{"inv_mult_quad"}, \code{"exponential"}, or \code{"poly"}.
#' @param ... Optional kernel parameters. Supported names are \code{gamma} and
#'   \code{poly_s}.
#' @return An \code{n x n} kernel matrix.
#' @export
compute_kernel_matrix_gpu = function(X, kernel = "gaussian", ...){
	args = list(...)
	device = if (!is.null(args$device)) args$device else 0
	if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))) {
		device = -1
	}
	gamma = if (!is.null(args$gamma)) args$gamma else 1
	poly_s = if (!is.null(args$poly_s)) args$poly_s else 2
	assertMatrix(X, mode = "numeric", any.missing = FALSE)
	assertChoice(kernel, c("gaussian", "laplacian", "inv_mult_quad", "exponential", "poly"))
	assertInt(device, lower = -1)
	assertCount(poly_s, positive = TRUE)
	assertNumeric(gamma, len = 1, lower = .Machine$double.eps, any.missing = FALSE)
	compute_kernel_matrix_gpu_cpp(X, kernel, as.integer(poly_s), gamma, as.integer(device))
}

#' Compute kernel objective values through the optional native backend
#'
#' @param W An \code{r x n} numeric matrix of design vectors.
#' @param Kgram An \code{n x n} kernel Gram matrix.
#' @param device Zero-based device id. -1 for native C++.
#' @return A numeric vector with one quadratic form per row of \code{W}.
#' @export
compute_objective_vals_gpu = function(W, Kgram, device = 0){
    if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))) {
        device = -1
    }
	assertMatrix(W, mode = "numeric", any.missing = FALSE)
	assertMatrix(Kgram, mode = "numeric", any.missing = FALSE)
	assertInt(device, lower = -1)
	compute_objective_vals_gpu_cpp(W, Kgram, as.integer(device))
}

#' Compute multiple kernel objective values using GPU
#'
#' @param W The design matrix (r x n)
#' @param Kgrams A list of Gram matrices
#' @param weights A vector of weights for the kernels
#' @param initial_objs A vector of initial objective values
#' @param running_sums A vector of current kernel sums
#' @param max_reds A vector of max reduction log objective values
#' @param maximum_gain_scaling Scaling factor
#' @param device The device ID
#' @export
compute_multiple_kernel_objective_vals_gpu = function(W, Kgrams, weights, initial_objs, running_sums, max_reds, maximum_gain_scaling, device = 0){
    if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))) {
        device = -1
    }
    compute_multiple_kernel_objective_vals_gpu_cpp(W, Kgrams, weights, initial_objs, running_sums, max_reds, maximum_gain_scaling, as.integer(device))
}


#' Compute randomization metrics through the optional native backend
#'
#' @param W An \code{r x n} integer matrix of design vectors (1/0).
#' @param device Zero-based device id. -1 for native C++.
#' @return An \code{n x n} matrix of same-group probabilities.
#' @export
compute_randomization_metrics_gpu = function(W, device = 0){
    if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))) {
        device = -1
    }
	assertMatrix(W, mode = "integer", any.missing = FALSE)
	assertInt(device, lower = -1)
	compute_randomization_metrics_gpu_cpp(W, as.integer(device))
}

#' Run a full greedy search on GPU (Upload Once, In-Place)
#'
#' @param X The design matrix
#' @param Sinv Inverse covariance matrix
#' @param start_indicT Starting allocation
#' @param max_iters Maximum iterations
#' @param device Device ID
#' @export
full_greedy_search_gpu = function(X, Sinv, start_indicT, max_iters = 100, device = 0) {
    if (!isTRUE(getOption("GreedyExperimentalDesign.use_gpu", FALSE))) {
        device = -1
    }
    full_greedy_search_gpu_cpp(X, Sinv, start_indicT, as.integer(max_iters), as.integer(device))
}
