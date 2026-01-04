#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool is_symmetric_cpp(NumericMatrix D, double tol) {
	int n = D.nrow();
	if (D.ncol() != n) {
		stop("D must be square");
	}
	if (tol < 0.0) {
		stop("tol must be non-negative");
	}
	for (int i = 0; i < n - 1; i++){
		for (int j = i + 1; j < n; j++){
			double a = D(i, j);
			double b = D(j, i);
			bool a_na = NumericVector::is_na(a);
			bool b_na = NumericVector::is_na(b);
			if (a_na && b_na) {
				continue;
			}
			if (a_na != b_na) {
				return false;
			}
			if (std::fabs(a - b) > tol) {
				return false;
			}
		}
	}
	return true;
}

struct Pair {
	int a;
	int b;
	bool is_na;
};

static bool pair_less(const Pair& lhs, const Pair& rhs) {
	if (lhs.is_na != rhs.is_na) {
		return rhs.is_na;
	}
	if (lhs.is_na) {
		return false;
	}
	if (lhs.a != rhs.a) {
		return lhs.a < rhs.a;
	}
	return lhs.b < rhs.b;
}

// [[Rcpp::export]]
IntegerMatrix sort_unique_pairs_cpp(IntegerMatrix pairs) {
	if (pairs.ncol() != 2) {
		stop("pairs must have 2 columns");
	}
	int n = pairs.nrow();
	std::vector<Pair> vec;
	vec.reserve(n);
	for (int i = 0; i < n; i++){
		int a = pairs(i, 0);
		int b = pairs(i, 1);
		if (a == NA_INTEGER || b == NA_INTEGER) {
			vec.push_back({NA_INTEGER, NA_INTEGER, true});
			continue;
		}
		if (a > b) {
			int tmp = a;
			a = b;
			b = tmp;
		}
		vec.push_back({a, b, false});
	}
	std::sort(vec.begin(), vec.end(), pair_less);
	std::vector<Pair> uniq;
	uniq.reserve(vec.size());
	for (size_t i = 0; i < vec.size(); i++){
		if (i == 0) {
			uniq.push_back(vec[i]);
			continue;
		}
		const Pair& prev = uniq.back();
		const Pair& cur = vec[i];
		bool same = false;
		if (prev.is_na && cur.is_na) {
			same = true;
		} else if (!prev.is_na && !cur.is_na && prev.a == cur.a && prev.b == cur.b) {
			same = true;
		}
		if (!same) {
			uniq.push_back(cur);
		}
	}

	IntegerMatrix out(uniq.size(), 2);
	for (size_t i = 0; i < uniq.size(); i++){
		if (uniq[i].is_na) {
			out(i, 0) = NA_INTEGER;
			out(i, 1) = NA_INTEGER;
		} else {
			out(i, 0) = uniq[i].a;
			out(i, 1) = uniq[i].b;
		}
	}
	return out;
}
