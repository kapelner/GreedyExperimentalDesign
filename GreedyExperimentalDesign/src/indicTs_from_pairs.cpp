#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix compute_indicTs_from_pairs_cpp(IntegerMatrix pairs, IntegerMatrix w_diff, int n) {
	int num_pairs = pairs.nrow();
	if (pairs.ncol() != 2) {
		stop("pairs must have 2 columns");
	}
	if (w_diff.ncol() != num_pairs) {
		stop("w_diff must have nrow(pairs) columns");
	}
	int r = w_diff.nrow();
	IntegerMatrix indicTs(r, n);
	for (int row = 0; row < r; row++){
		for (int i = 0; i < num_pairs; i++){
			int a = pairs(i, 0) - 1;
			int b = pairs(i, 1) - 1;
			if (a < 0 || b < 0 || a >= n || b >= n) {
				stop("pairs indices out of bounds");
			}
			if (w_diff(row, i) == 1) {
				indicTs(row, a) = 1;
				indicTs(row, b) = 0;
			} else {
				indicTs(row, a) = 0;
				indicTs(row, b) = 1;
			}
		}
	}
	return indicTs;
}

// [[Rcpp::export]]
IntegerMatrix compute_indicTs_from_pairs_lgl_cpp(IntegerMatrix pairs, LogicalMatrix w_diff, int n) {
	int num_pairs = pairs.nrow();
	if (pairs.ncol() != 2) {
		stop("pairs must have 2 columns");
	}
	if (w_diff.ncol() != num_pairs) {
		stop("w_diff must have nrow(pairs) columns");
	}
	int r = w_diff.nrow();
	IntegerMatrix indicTs(r, n);
	for (int row = 0; row < r; row++){
		for (int i = 0; i < num_pairs; i++){
			int a = pairs(i, 0) - 1;
			int b = pairs(i, 1) - 1;
			if (a < 0 || b < 0 || a >= n || b >= n) {
				stop("pairs indices out of bounds");
			}
			int val = w_diff(row, i);
			if (val == NA_LOGICAL) {
				stop("w_diff contains NA");
			}
			if (val) {
				indicTs(row, a) = 1;
				indicTs(row, b) = 0;
			} else {
				indicTs(row, a) = 0;
				indicTs(row, b) = 1;
			}
		}
	}
	return indicTs;
}

// [[Rcpp::export]]
IntegerMatrix compute_indicTs_from_pairs_raw_cpp(IntegerMatrix pairs, RawMatrix w_diff, int n) {
	int num_pairs = pairs.nrow();
	if (pairs.ncol() != 2) {
		stop("pairs must have 2 columns");
	}
	if (w_diff.ncol() != num_pairs) {
		stop("w_diff must have nrow(pairs) columns");
	}
	int r = w_diff.nrow();
	IntegerMatrix indicTs(r, n);
	for (int row = 0; row < r; row++){
		for (int i = 0; i < num_pairs; i++){
			int a = pairs(i, 0) - 1;
			int b = pairs(i, 1) - 1;
			if (a < 0 || b < 0 || a >= n || b >= n) {
				stop("pairs indices out of bounds");
			}
			Rbyte val = w_diff(row, i);
			if (val == 1) {
				indicTs(row, a) = 1;
				indicTs(row, b) = 0;
			} else if (val == 0) {
				indicTs(row, a) = 0;
				indicTs(row, b) = 1;
			} else {
				stop("w_diff raw values must be 0 or 1");
			}
		}
	}
	return indicTs;
}
