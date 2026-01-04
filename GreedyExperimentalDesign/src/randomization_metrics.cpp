#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_randomization_metrics_cpp(NumericMatrix designs) {
	int n = designs.nrow();
	int r = designs.ncol();
	NumericMatrix p_hat_same_group(n, n);
	for (int i1 = 0; i1 < n - 1; i1++){
		for (int i2 = i1 + 1; i2 < n; i2++){
			int num_same_group = 0;
			for (int j = 0; j < r; j++){
				num_same_group += (designs(i1, j) == designs(i2, j)) ? 1 : 0;
			}
			p_hat_same_group(i1, i2) = num_same_group / static_cast<double>(r);
		}
	}

	double s_n = (n - 2) / (static_cast<double>(2 * n - 2));
	int num_pairs = n * (n - 1) / 2;

	double sum_entropies = 0.0;
	double sum_sqd_dev = 0.0;
	for (int i1 = 0; i1 < n - 1; i1++){
		for (int i2 = i1 + 1; i2 < n; i2++){
			double p_hat = p_hat_same_group(i1, i2);
			if (p_hat > 0.0) {
				sum_entropies += p_hat * std::log(p_hat);
			}
			double one_minus = 1.0 - p_hat;
			if (one_minus > 0.0) {
				sum_entropies += one_minus * std::log(one_minus);
			}
			double diff = p_hat - s_n;
			sum_sqd_dev += diff * diff;
		}
	}

	double entropy_norm_factor = s_n * std::log(s_n) + (1.0 - s_n) * std::log(1.0 - s_n);
	double entropy_metric = (1.0 / static_cast<double>(num_pairs)) * (sum_entropies / entropy_norm_factor);
	double const_factor = (2.0 / static_cast<double>(n)) * std::sqrt((2.0 * n - 2.0) / static_cast<double>(n - 2));
	double se_metric = const_factor * std::sqrt(sum_sqd_dev);

	return List::create(
		Named("p_hat_ijs") = p_hat_same_group,
		Named("rand_entropy_metric") = entropy_metric,
		Named("rand_norm_se_metric") = se_metric
	);
}
