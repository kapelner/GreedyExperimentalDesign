#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// kernel_type: 1=poly, 2=exponential, 3=laplacian, 4=inv_mult_quad, 5=gaussian
// [[Rcpp::export]]
NumericMatrix compute_kernel_matrix_cpp(NumericMatrix X, int kernel_type, int poly_s) {
	int n = X.nrow();
	int p = X.ncol();
	if (kernel_type < 1 || kernel_type > 5) {
		stop("Invalid kernel_type");
	}
	if (kernel_type == 1 && poly_s <= 0) {
		stop("poly_s must be positive for polynomial kernel");
	}
	NumericMatrix K(n, n);
	for (int i = 0; i < n; i++){
		if (kernel_type == 1 || kernel_type == 2) {
			double dot = 0.0;
			for (int k = 0; k < p; k++){
				dot += X(i, k) * X(i, k);
			}
			K(i, i) = (kernel_type == 1) ? std::pow(1.0 + dot / poly_s, poly_s) : std::exp(dot);
		} else {
			K(i, i) = 1.0;
		}
		for (int j = i + 1; j < n; j++){
			double val = 0.0;
			if (kernel_type == 1 || kernel_type == 2){
				double dot = 0.0;
				for (int k = 0; k < p; k++){
					dot += X(i, k) * X(j, k);
				}
				if (kernel_type == 1){
					val = std::pow(1.0 + dot / poly_s, poly_s);
				} else {
					val = std::exp(dot);
				}
			} else {
				double sqd = 0.0;
				for (int k = 0; k < p; k++){
					double diff = X(i, k) - X(j, k);
					sqd += diff * diff;
				}
				if (kernel_type == 3){
					val = std::exp(-std::sqrt(sqd));
				} else if (kernel_type == 4){
					val = 1.0 / std::sqrt(sqd + 1.0);
				} else {
					val = std::exp(-sqd);
				}
			}
			K(i, j) = val;
			K(j, i) = val;
		}
	}
	return K;
}
