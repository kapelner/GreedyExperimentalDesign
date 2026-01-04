#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List greedy_orthogonalization_eliminate_cpp(NumericMatrix Rij, int Rmin, bool stop_if_worse) {
	int R = Rij.nrow();
	if (R != Rij.ncol()) {
		stop("Rij must be square");
	}
	if (Rmin < 2 || Rmin > R) {
		stop("Rmin must be between 2 and nrow(Rij)");
	}

	std::vector<double> row_sum(R, 0.0);
	for (int i = 0; i < R; i++){
		double sum = 0.0;
		for (int j = 0; j < R; j++){
			if (i == j) {
				continue;
			}
			double val = Rij(i, j);
			if (!NumericVector::is_na(val)) {
				sum += val;
			}
		}
		row_sum[i] = sum;
	}

	double total_offdiag = 0.0;
	for (int i = 0; i < R; i++){
		total_offdiag += row_sum[i];
	}

	int active_count = R;
	std::vector<int> eliminated;
	eliminated.reserve(R - Rmin);
	NumericVector avg_abs_rijss;
	avg_abs_rijss.push_back(total_offdiag / (static_cast<double>(R) * (R - 1)));

	std::vector<char> active(R, 1);

	while (active_count > Rmin) {
		double max_sum = R_NegInf;
		int left_star = -1;
		for (int i = 0; i < R; i++){
			if (!active[i]) {
				continue;
			}
			if (row_sum[i] > max_sum) {
				max_sum = row_sum[i];
				left_star = i;
			}
		}
		if (left_star < 0) {
			break;
		}

		double new_total = total_offdiag - 2.0 * row_sum[left_star];
		double denom = static_cast<double>(active_count - 1) * (active_count - 2);
		double new_avg = new_total / denom;
		if (stop_if_worse && new_avg > avg_abs_rijss[avg_abs_rijss.size() - 1]) {
			break;
		}

		eliminated.push_back(left_star + 1);
		avg_abs_rijss.push_back(new_avg);
		total_offdiag = new_total;

		active[left_star] = 0;
		for (int j = 0; j < R; j++){
			if (!active[j]) {
				continue;
			}
			double val = Rij(left_star, j);
			if (!NumericVector::is_na(val)) {
				row_sum[j] -= val;
			}
		}
		row_sum[left_star] = NA_REAL;
		active_count--;
	}

	IntegerVector remaining(active_count);
	int idx = 0;
	for (int i = 0; i < R; i++){
		if (active[i]) {
			remaining[idx++] = i + 1;
		}
	}

	IntegerVector eliminated_vec(eliminated.begin(), eliminated.end());
	return List::create(
		Named("remaining") = remaining,
		Named("eliminated") = eliminated_vec,
		Named("avg_abs_rijss") = avg_abs_rijss
	);
}
