#include <Rcpp.h>
#include <random>
#include <chrono>
using namespace Rcpp;

static std::default_random_engine init_rng(int seed) {
	unsigned int rng_seed;
	if (seed == NA_INTEGER) {
		rng_seed = static_cast<unsigned int>(
			std::chrono::system_clock::now().time_since_epoch().count()
		);
	} else {
		rng_seed = static_cast<unsigned int>(seed);
	}
	return std::default_random_engine(rng_seed);
}

static void shuffle_in_place(std::vector<int>& vec, std::default_random_engine& rng) {
	for (int i = static_cast<int>(vec.size()) - 1; i > 0; i--){
		std::uniform_int_distribution<int> dist(0, i);
		int j = dist(rng);
		std::swap(vec[i], vec[j]);
	}
}

// [[Rcpp::export]]
IntegerMatrix complete_randomization_forced_balanced_cpp(int n, int r, int seed) {
	std::default_random_engine rng = init_rng(seed);
	IntegerMatrix indicTs(r, n);
	std::vector<int> vec(n, 0);
	int nT = n / 2;
	for (int i = 0; i < nT; i++){
		vec[i] = 1;
	}
	for (int row = 0; row < r; row++){
		shuffle_in_place(vec, rng);
		for (int j = 0; j < n; j++){
			indicTs(row, j) = vec[j];
		}
	}
	return indicTs;
}

// [[Rcpp::export]]
IntegerMatrix complete_randomization_imbalanced_cpp(int n, int nT, int r, int seed) {
	std::default_random_engine rng = init_rng(seed);
	IntegerMatrix indicTs(r, n);
	std::vector<int> vec(n, 0);
	for (int i = 0; i < nT; i++){
		vec[i] = 1;
	}
	for (int row = 0; row < r; row++){
		shuffle_in_place(vec, rng);
		for (int j = 0; j < n; j++){
			indicTs(row, j) = vec[j];
		}
	}
	return indicTs;
}

// [[Rcpp::export]]
IntegerMatrix complete_randomization_cpp(int n, int r) {
	IntegerMatrix indicTs(r, n);
	for (int row = 0; row < r; row++){
		for (int j = 0; j < n; j++){
			indicTs(row, j) = (R::runif(0.0, 1.0) < 0.5) ? 1 : 0;
		}
	}
	return indicTs;
}

// [[Rcpp::export]]
IntegerMatrix imbalanced_block_designs_cpp(int n_B, int n_B_T, int B, int r, int seed) {
	std::default_random_engine rng = init_rng(seed);
	IntegerMatrix indicTs(r, n_B * B);
	std::vector<int> vec(n_B, 0);
	for (int i = 0; i < n_B_T; i++){
		vec[i] = 1;
	}
	for (int b = 0; b < B; b++){
		int col_offset = b * n_B;
		for (int row = 0; row < r; row++){
			shuffle_in_place(vec, rng);
			for (int j = 0; j < n_B; j++){
				indicTs(row, col_offset + j) = vec[j];
			}
		}
	}
	return indicTs;
}
