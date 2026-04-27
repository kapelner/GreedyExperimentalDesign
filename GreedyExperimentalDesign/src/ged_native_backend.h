#ifndef GED_NATIVE_BACKEND_H
#define GED_NATIVE_BACKEND_H

#include <Rcpp.h>

namespace ged_native_backend {

Rcpp::LogicalVector gpu_available();
Rcpp::DataFrame gpu_devices();
Rcpp::NumericMatrix distance_matrix(Rcpp::NumericMatrix X,
                                   const std::string& backend,
                                   int device);
Rcpp::NumericMatrix kernel_matrix(Rcpp::NumericMatrix X,
                                  const std::string& kernel,
                                  int poly_s,
                                  double gamma,
                                  int device);
Rcpp::NumericVector objective_vals(Rcpp::NumericMatrix W,
		Rcpp::NumericMatrix Kgram, int device);
Rcpp::NumericVector multiple_kernel_objective_vals(Rcpp::NumericMatrix W,
        Rcpp::List Kgrams, Rcpp::NumericVector weights,
        Rcpp::NumericVector initial_objs, Rcpp::NumericVector running_sums,
        Rcpp::NumericVector max_reds, double maximum_gain_scaling, int device);
Rcpp::IntegerVector full_greedy_search(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Sinv, Rcpp::IntegerVector start_indicT, int max_iters, int device);
Rcpp::NumericMatrix randomization_metrics(Rcpp::IntegerMatrix W, int device);


}

#endif
