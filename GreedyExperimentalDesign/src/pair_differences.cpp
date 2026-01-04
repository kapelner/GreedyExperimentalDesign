#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_pair_differences_cpp(NumericMatrix X, IntegerMatrix pairs) {
	int num_pairs = pairs.nrow();
	if (pairs.ncol() != 2) {
		stop("pairs must have 2 columns");
	}
	int p = X.ncol();
	NumericMatrix diffs(num_pairs, p);
	for (int i = 0; i < num_pairs; i++){
		int a = pairs(i, 0) - 1;
		int b = pairs(i, 1) - 1;
		if (a < 0 || b < 0 || a >= X.nrow() || b >= X.nrow()) {
			stop("pairs indices out of bounds");
		}
		for (int j = 0; j < p; j++){
			diffs(i, j) = X(a, j) - X(b, j);
		}
	}
	return diffs;
}
