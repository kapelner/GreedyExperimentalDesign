#include <Rcpp.h>
#include <random>
#include <chrono>
using namespace Rcpp;

// [[Rcpp::export]]
List generate_block_design_cpp(int B, int nR, NumericVector dummy_block) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	List Ws(B);
	int n_B = dummy_block.size();
	for (int b = 0; b < B; b++){
	  NumericMatrix Ws_b(nR, n_B);
	  for (int nr = 0; nr < nR; nr++){
	    std::shuffle(dummy_block.begin(), dummy_block.end(), generator);
	    Ws_b(nr, _) = dummy_block;
	  }
	  Ws[b] = Ws_b;
	}
	return Ws;
}
