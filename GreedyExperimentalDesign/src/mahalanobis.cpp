#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_mahalanobis_distance_matrix_cpp(NumericMatrix X, NumericMatrix Sinv) {
	int n = X.nrow();
	int p = X.ncol();
	if (Sinv.nrow() != p || Sinv.ncol() != p) {
		stop("Sinv must be a square matrix with dimensions matching ncol(X)");
	}
	NumericMatrix D(n, n);
	std::fill(D.begin(), D.end(), NA_REAL);
	std::vector<double> diff(p);
	for (int i_1 = 0; i_1 < (n - 1); i_1++){
		for (int i_2 = i_1 + 1; i_2 < n; i_2++){
			for (int j = 0; j < p; j++){
				diff[j] = X(i_1, j) - X(i_2, j);
			}
			double dist = 0.0;
			for (int j = 0; j < p; j++){
				double inner = 0.0;
				for (int k = 0; k < p; k++){
					inner += Sinv(j, k) * diff[k];
				}
				dist += diff[j] * inner;
			}
			D(i_1, i_2) = dist;
			D(i_2, i_1) = dist;
		}
	}
	return D;
}
