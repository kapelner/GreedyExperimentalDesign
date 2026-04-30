#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix standardize_data_matrix_cpp(NumericMatrix X) {
	int n = X.nrow();
	int p = X.ncol();
	NumericMatrix Xstd(n, p);
	for (int j = 0; j < p; j++){
		bool has_na = false;
		double sum = 0.0;
		double sumsq = 0.0;
		for (int i = 0; i < n; i++){
			double x = X(i, j);
			if (NumericVector::is_na(x)) {
				has_na = true;
				break;
			}
			sum += x;
			sumsq += x * x;
		}
		if (has_na || n <= 1) {
			for (int i = 0; i < n; i++){
				Xstd(i, j) = NA_REAL;
			}
			continue;
		}
		double mean = sum / static_cast<double>(n);
		double var = (sumsq - (sum * sum) / static_cast<double>(n)) / static_cast<double>(n - 1);
		if (var < 0.0) {
			var = 0.0;
		}
		double sd = std::sqrt(var);
		for (int i = 0; i < n; i++){
			Xstd(i, j) = (X(i, j) - mean) / sd;
		}
	}
	return Xstd;
}

// [[Rcpp::export]]
double compute_objective_val_cpp(NumericMatrix X, IntegerVector indic_T, std::string objective, Nullable<NumericMatrix> inv_cov_X = R_NilValue) {
	int n = X.nrow();
	int p = X.ncol();
	if (indic_T.size() != n) {
		stop("indic_T length must match nrow(X)");
	}
	std::vector<double> sum_T(p, 0.0);
	std::vector<double> sum_C(p, 0.0);
	std::vector<double> sum_all(p, 0.0);
	std::vector<double> sumsq_all(p, 0.0);
	int nT = 0;
	int nC = 0;

	for (int i = 0; i < n; i++){
		int t = indic_T[i];
		if (t != 0 && t != 1) {
			stop("indic_T must be binary");
		}
		if (t == 1) {
			nT++;
		} else {
			nC++;
		}
		for (int j = 0; j < p; j++){
			double x = X(i, j);
			sum_all[j] += x;
			sumsq_all[j] += x * x;
			if (t == 1) {
				sum_T[j] += x;
			} else {
				sum_C[j] += x;
			}
		}
	}
	if (nT == 0 || nC == 0) {
		stop("Both treatment and control groups must be non-empty");
	}

	std::vector<double> diff(p, 0.0);
	for (int j = 0; j < p; j++){
		diff[j] = (sum_T[j] / static_cast<double>(nT)) - (sum_C[j] / static_cast<double>(nC));
	}

	if (objective == "abs_sum_diff") {
		if (n <= 1) {
			return NA_REAL;
		}
		double total = 0.0;
		for (int j = 0; j < p; j++){
			double var = (sumsq_all[j] - (sum_all[j] * sum_all[j]) / static_cast<double>(n)) / static_cast<double>(n - 1);
			if (var < 0.0) {
				var = 0.0;
			}
			double sd = std::sqrt(var);
			total += std::fabs(diff[j] / sd);
		}
		return total;
	}

	if (objective == "mahal_dist") {
		if (inv_cov_X.isNull()) {
			stop("inv_cov_X is required for mahal_dist");
		}
		NumericMatrix Sinv(inv_cov_X);
		if (Sinv.nrow() != p || Sinv.ncol() != p) {
			stop("inv_cov_X must be p x p");
		}
		double val = 0.0;
		for (int i = 0; i < p; i++){
			for (int k = 0; k < p; k++){
				val += diff[i] * Sinv(i, k) * diff[k];
			}
		}
		return val;
	}

	stop("objective invalid");
	return NA_REAL;
}
