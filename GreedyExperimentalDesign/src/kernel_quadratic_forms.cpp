#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_kernel_quadratic_forms_cpp(NumericMatrix W, NumericMatrix K) {
	int r = W.nrow();
	int n = W.ncol();
	if (K.nrow() != n || K.ncol() != n) {
		stop("K must be n x n with n = ncol(W)");
	}
	NumericVector vals(r);
	for (int row = 0; row < r; row++){
		double sum = 0.0;
		for (int i = 0; i < n; i++){
			double wi = W(row, i);
			for (int j = 0; j < n; j++){
				sum += wi * K(i, j) * W(row, j);
			}
		}
		vals[row] = sum;
	}
	return vals;
}
