#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gen_var_cov_matrix_block_designs_cpp(int n, double prop_T, int B) {
	if (B <= 0) {
		stop("B must be positive");
	}
	if (n <= 0) {
		stop("n must be positive");
	}
	int n_B = n / B;
	if (n_B * B != n) {
		stop("n must be divisible by B");
	}
	if (n_B <= 1) {
		stop("n/B must be greater than 1");
	}
	double n_T = n * prop_T;
	double n_C = n * (1.0 - prop_T);
	double b = (n_T - n_C) * (n_T - n_C) / (static_cast<double>(n) * n);
	double diag_val = 1.0 - b;
	double off_val = -(1.0 - b) / static_cast<double>(n_B - 1);

	NumericMatrix SigmaW(n, n);
	for (int block = 0; block < B; block++){
		int start = block * n_B;
		for (int i = 0; i < n_B; i++){
			for (int j = 0; j < n_B; j++){
				SigmaW(start + i, start + j) = (i == j) ? diag_val : off_val;
			}
		}
	}
	return SigmaW;
}
