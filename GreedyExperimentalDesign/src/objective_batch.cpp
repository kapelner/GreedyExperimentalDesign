#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_objective_vals_cpp(NumericMatrix X, IntegerMatrix indicTs, std::string objective, Nullable<NumericMatrix> inv_cov_X = R_NilValue) {
	int n = X.nrow();
	int p = X.ncol();
	int r = indicTs.nrow();
	if (indicTs.ncol() != n) {
		stop("indicTs must have n columns matching nrow(X)");
	}

	std::vector<double> sum_all(p, 0.0);
	std::vector<double> sumsq_all(p, 0.0);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < p; j++){
			double x = X(i, j);
			sum_all[j] += x;
			sumsq_all[j] += x * x;
		}
	}

	std::vector<double> sd_all;
	if (objective == "abs_sum_diff") {
		if (n <= 1) {
			stop("n must be greater than 1 for abs_sum_diff");
		}
		sd_all.resize(p);
		for (int j = 0; j < p; j++){
			double var = (sumsq_all[j] - (sum_all[j] * sum_all[j]) / static_cast<double>(n)) / static_cast<double>(n - 1);
			if (var < 0.0) {
				var = 0.0;
			}
			sd_all[j] = std::sqrt(var);
		}
	} else if (objective == "mahal_dist") {
		if (inv_cov_X.isNull()) {
			stop("inv_cov_X is required for mahal_dist");
		}
	} else {
		stop("objective invalid");
	}

	NumericMatrix Sinv;
	if (objective == "mahal_dist") {
		Sinv = NumericMatrix(inv_cov_X);
		if (Sinv.nrow() != p || Sinv.ncol() != p) {
			stop("inv_cov_X must be p x p");
		}
	}

	NumericVector vals(r);
	std::vector<double> sum_T(p);
	std::vector<double> diff(p);

	for (int row = 0; row < r; row++){
		std::fill(sum_T.begin(), sum_T.end(), 0.0);
		int nT = 0;
		for (int i = 0; i < n; i++){
			int t = indicTs(row, i);
			if (t != 0 && t != 1) {
				stop("indicTs must be binary");
			}
			if (t == 1) {
				nT++;
				for (int j = 0; j < p; j++){
					sum_T[j] += X(i, j);
				}
			}
		}
		int nC = n - nT;
		if (nT == 0 || nC == 0) {
			stop("Both treatment and control groups must be non-empty");
		}

		for (int j = 0; j < p; j++){
			double mean_T = sum_T[j] / static_cast<double>(nT);
			double mean_C = (sum_all[j] - sum_T[j]) / static_cast<double>(nC);
			diff[j] = mean_T - mean_C;
		}

		if (objective == "abs_sum_diff") {
			double total = 0.0;
			for (int j = 0; j < p; j++){
				total += std::fabs(diff[j] / sd_all[j]);
			}
			vals[row] = total;
		} else {
			double val = 0.0;
			for (int i = 0; i < p; i++){
				for (int k = 0; k < p; k++){
					val += diff[i] * Sinv(i, k) * diff[k];
				}
			}
			vals[row] = val;
		}
	}

	return vals;
}
