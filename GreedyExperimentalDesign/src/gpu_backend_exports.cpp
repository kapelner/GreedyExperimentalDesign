#include <Rcpp.h>
#include <string>
#include "ged_native_backend.h"

#ifdef GED_USE_WGPU
namespace ged_wgpu_backend { void release_context(); }
#endif

using namespace Rcpp;

// [[Rcpp::export]]
void ged_gpu_release_context_cpp() {
#ifdef GED_USE_WGPU
    ged_wgpu_backend::release_context();
#endif
}

// [[Rcpp::export]]
LogicalVector ged_gpu_available_cpp() {
	return ged_native_backend::gpu_available();
}

// [[Rcpp::export]]
DataFrame ged_gpu_devices_cpp() {
	return ged_native_backend::gpu_devices();
}

// [[Rcpp::export]]
NumericMatrix compute_distance_matrix_gpu_cpp(NumericMatrix X,
		std::string backend = "auto", int device = 0) {
	return ged_native_backend::distance_matrix(X, backend, device);
}

// [[Rcpp::export]]
NumericMatrix compute_kernel_matrix_gpu_cpp(NumericMatrix X,
		std::string kernel = "gaussian", int poly_s = 2,
		double gamma = 1.0, int device = 0) {
	return ged_native_backend::kernel_matrix(X, kernel, poly_s, gamma, device);
}

// [[Rcpp::export]]
NumericVector compute_objective_vals_gpu_cpp(NumericMatrix W,
		NumericMatrix Kgram, int device = 0) {
	return ged_native_backend::objective_vals(W, Kgram, device);
}

// [[Rcpp::export]]
NumericVector compute_multiple_kernel_objective_vals_gpu_cpp(NumericMatrix W,
		List Kgrams, NumericVector weights, NumericVector initial_objs,
		NumericVector running_sums, NumericVector max_reds,
		double maximum_gain_scaling, int device = 0) {
	return ged_native_backend::multiple_kernel_objective_vals(W, Kgrams, weights, initial_objs, running_sums, max_reds, maximum_gain_scaling, device);
}

// [[Rcpp::export]]
NumericMatrix compute_randomization_metrics_gpu_cpp(IntegerMatrix W, int device = 0) {
	return ged_native_backend::randomization_metrics(W, device);
}

// [[Rcpp::export]]
IntegerVector full_greedy_search_gpu_cpp(NumericMatrix X, NumericMatrix Sinv, IntegerVector start_indicT, int max_iters, int device = 0) {
	return ged_native_backend::full_greedy_search(X, Sinv, start_indicT, max_iters, device);
}
