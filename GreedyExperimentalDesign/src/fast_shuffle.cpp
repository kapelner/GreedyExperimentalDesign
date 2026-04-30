#include <Rcpp.h>
#include <random>
#include <chrono>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector shuffle_cpp(NumericVector w, unsigned int seed) {
	//check if seed is null; if so, set it below from the system clock
	if ((int)seed == NA_INTEGER){
		seed = std::chrono::system_clock::now().time_since_epoch().count();
	}
	NumericVector w_shuffled = clone(w);
	std::shuffle(w_shuffled.begin(), w_shuffled.end(), std::default_random_engine(seed));
	return w_shuffled;
}

// [[Rcpp::export]]
LogicalMatrix create_all_ys_cpp(NumericVector pCs, NumericVector pTs, NumericMatrix W, int two_n, int nY){
  LogicalMatrix Y(nY, two_n);
  for (int i = 0; i < nY; i++){
    for (int j = 0; j < two_n; j++){
      double u = R::runif(0.0, 1.0);
      if (W(i, j) == 1){
        Y(i, j) = (u <= pTs(j));
      } else {
        Y(i, j) = (u <= pCs(j));
      }
    }
  }
  return Y;
} 

// [[Rcpp::export]]
bool all_elements_same_cpp(NumericVector vec) {
	if (vec.size() <= 1) {
	  return true;
	}
	double firstElement = vec[0];
	for (R_xlen_t i = 1; i < vec.size(); ++i) {
	  if (vec[i] != firstElement) {
	    return false;
	  }
	}
	return true;
}  

// [[Rcpp::export]]
NumericMatrix gen_pm_designs_cpp(NumericMatrix indicies_pairs, int n, int r){
	NumericMatrix W(r, n * 2);
	for (int w = 0; w < r; w++){
	  for (int i = 0; i < n; i++){
	    int a = (R::runif(0.0, 1.0) < 0.5) ? -1 : 1;
	    int idx0 = static_cast<int>(indicies_pairs(i, 0)) - 1;
	    int idx1 = static_cast<int>(indicies_pairs(i, 1)) - 1;
	    if (idx0 < 0 || idx1 < 0 || idx0 >= n * 2 || idx1 >= n * 2){
	      stop("indicies_pairs contains out-of-range indices.");
	    }
	    W(w, idx0) = a;
	    W(w, idx1) = -a;
	  }
	}
	return W;
}
