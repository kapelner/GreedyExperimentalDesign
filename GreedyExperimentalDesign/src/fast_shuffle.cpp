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
	std::shuffle(w.begin(), w.end(), std::default_random_engine(seed));
	return w;
}

// [[Rcpp::export]]
LogicalMatrix create_all_ys_cpp(NumericVector pCs, NumericVector pTs, NumericMatrix W, int two_n, int nY){
  LogicalMatrix Y(nY, two_n);
  for (int i = 0; i < nY; i++){
    for (int j = 0; j < two_n; j++){
      if (W(i, j) == 1){
        Y(i, j) = (rand() / (RAND_MAX + 1.0) <= pTs(j));
      } else {
        Y(i, j) = (rand() / (RAND_MAX + 1.0) <= pCs(j));
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
	int firstElement = vec[0];
	for (size_t i = 1; i < vec.size(); ++i) {
	  if (vec[i] != firstElement) {
	    return false;
	  }
	}
	return true;
}  