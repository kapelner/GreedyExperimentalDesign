#include "ged_native_backend.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include "ged_gpu_config.h"

#ifdef GED_USE_WGPU
namespace ged_wgpu_backend {
    Rcpp::NumericMatrix distance_matrix(Rcpp::NumericMatrix X, int device_id);
    Rcpp::NumericMatrix kernel_matrix(Rcpp::NumericMatrix X, const std::string& kernel, int poly_s, double gamma, int device_id);
    Rcpp::NumericVector objective_vals(Rcpp::NumericMatrix W, Rcpp::NumericMatrix Kgram, int device_id);
    Rcpp::NumericVector multiple_kernel_objective_vals(Rcpp::NumericMatrix W,
            Rcpp::List Kgrams, Rcpp::NumericVector weights,
            Rcpp::NumericVector initial_objs, Rcpp::NumericVector running_sums,
            Rcpp::NumericVector max_reds, double maximum_gain_scaling, int device);
    Rcpp::IntegerVector full_greedy_search(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Sinv, Rcpp::IntegerVector start_indicT, int max_iters, int device_id);
}
#endif


namespace {

struct DeviceInfo {
	std::string name;
	std::string backend;
	std::string type;
	bool usable;
};

std::string normalize_backend(const std::string& backend) {
	std::string value = backend;
	std::transform(value.begin(), value.end(), value.begin(), ::tolower);
	if (value != "auto" && value != "cpu" && value != "gpu") {
		Rcpp::stop("backend must be one of 'auto', 'cpu', or 'gpu'");
	}
	return value;
}

std::string normalize_kernel(const std::string& kernel) {
	std::string value = kernel;
	std::transform(value.begin(), value.end(), value.begin(), ::tolower);
	if (value != "gaussian" && value != "laplacian" &&
			value != "inv_mult_quad" && value != "exponential" &&
			value != "poly") {
		Rcpp::stop("kernel must be one of 'gaussian', 'laplacian', 'inv_mult_quad', 'exponential', or 'poly'");
	}
	return value;
}

Rcpp::NumericMatrix distance_matrix_cpu(Rcpp::NumericMatrix X) {
	const int n = X.nrow();
	const int p = X.ncol();
	Rcpp::NumericMatrix D(n, n);
	std::fill(D.begin(), D.end(), NA_REAL);
	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			double sqd = 0.0;
			for (int k = 0; k < p; ++k) {
				const double diff = X(i, k) - X(j, k);
				sqd += diff * diff;
			}
			D(i, j) = sqd;
			D(j, i) = sqd;
		}
	}
	return D;
}

Rcpp::NumericMatrix kernel_matrix_cpu(Rcpp::NumericMatrix X,
		const std::string& kernel, int poly_s, double gamma) {
	const std::string kernel_name = normalize_kernel(kernel);
	const int n = X.nrow();
	const int p = X.ncol();
	Rcpp::NumericMatrix K(n, n);
	for (int i = 0; i < n; ++i) {
        if (kernel_name == "poly" || kernel_name == "exponential") {
            double dot = 0.0;
            for (int k = 0; k < p; k++) dot += X(i, k) * X(i, k);
            K(i, i) = (kernel_name == "poly") ? std::pow(1.0 + dot / poly_s, poly_s) : std::exp(gamma * dot);
        } else {
		    K(i, i) = 1.0;
        }
		for (int j = i + 1; j < n; ++j) {
			double val = 0.0;
            if (kernel_name == "poly" || kernel_name == "exponential") {
                double dot = 0.0;
                for (int k = 0; k < p; k++) dot += X(i, k) * X(j, k);
                val = (kernel_name == "poly") ? std::pow(1.0 + dot / poly_s, poly_s) : std::exp(gamma * dot);
            } else {
                double sqd = 0.0;
                for (int k = 0; k < p; ++k) {
                    const double diff = X(i, k) - X(j, k);
                    sqd += diff * diff;
                }
                if (kernel_name == "laplacian") {
                    val = std::exp(-gamma * std::sqrt(sqd));
                } else if (kernel_name == "inv_mult_quad") {
                    val = 1.0 / std::sqrt(gamma * sqd + 1.0);
                } else {
                    val = std::exp(-gamma * sqd);
                }
            }
			K(i, j) = val;
			K(j, i) = val;
		}
	}
	return K;
}

Rcpp::NumericVector objective_vals_cpu(Rcpp::NumericMatrix W,
		Rcpp::NumericMatrix K) {
	const int r = W.nrow();
	const int n = W.ncol();
	Rcpp::NumericVector vals(r);
	for (int row = 0; row < r; ++row) {
		double sum = 0.0;
		for (int i = 0; i < n; ++i) {
			const double wi = W(row, i);
			for (int j = 0; j < n; ++j) {
				sum += wi * K(i, j) * W(row, j);
			}
		}
		vals[row] = sum;
	}
	return vals;
}

Rcpp::NumericVector multiple_kernel_objective_vals_cpu(Rcpp::NumericMatrix W,
		Rcpp::List Kgrams, Rcpp::NumericVector weights,
		Rcpp::NumericVector initial_objs, Rcpp::NumericVector running_sums,
		Rcpp::NumericVector max_reds, double maximum_gain_scaling) {
    int r = W.nrow();
    int m = Kgrams.size();
    Rcpp::NumericVector final_vals(r);
    for (int k = 0; k < m; k++) {
        Rcpp::NumericMatrix K = Kgrams[k];
        Rcpp::NumericVector kv = objective_vals_cpu(W, K);
        double w = weights[k];
        for (int i = 0; i < r; i++) {
            double val = kv[i];
            if (maximum_gain_scaling) {
                val = std::log10(initial_objs[k] / val) / max_reds[k];
            }
            final_vals[i] += w * val;
        }
    }
    return final_vals;
}

Rcpp::NumericMatrix randomization_metrics_cpu(Rcpp::IntegerMatrix W) {
    const int r = W.nrow();
    const int n = W.ncol();
    Rcpp::NumericMatrix P(n, n);
    for (int i1 = 0; i1 < n; i1++) {
        P(i1, i1) = 1.0;
        for (int i2 = i1 + 1; i2 < n; i2++) {
            int count = 0;
            for (int k = 0; k < r; k++) {
                if (W(k, i1) == W(k, i2)) {
                    count++;
                }
            }
            double val = (double)count / r;
            P(i1, i2) = val;
            P(i2, i1) = val;
        }
    }
    return P;
}

}

namespace ged_native_backend {

Rcpp::LogicalVector gpu_available() {
    bool has_gpu = false;
#if defined(GED_USE_WGPU)
    has_gpu = true;
#endif
	Rcpp::LogicalVector out = Rcpp::LogicalVector::create(has_gpu);
	out.attr("ggml_compiled") = false;
#if defined(GED_USE_WGPU)
	out.attr("wgpu_compiled") = true;
#else
	out.attr("wgpu_compiled") = false;
#endif
	return out;
}

Rcpp::DataFrame gpu_devices() {
	std::vector<DeviceInfo> devices;
#if defined(GED_USE_WGPU)
    devices.push_back(DeviceInfo{"WebGPU Device", "webgpu", "gpu", true});
#endif
	Rcpp::IntegerVector id(devices.size());
	Rcpp::CharacterVector name(devices.size());
	Rcpp::CharacterVector backend(devices.size());
	Rcpp::CharacterVector type(devices.size());
	Rcpp::LogicalVector usable(devices.size());
	for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(devices.size()); ++i) {
		id[i] = static_cast<int>(i);
		name[i] = devices[i].name;
		backend[i] = devices[i].backend;
		type[i] = devices[i].type;
		usable[i] = devices[i].usable;
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("id") = id,
		Rcpp::Named("name") = name,
		Rcpp::Named("backend") = backend,
		Rcpp::Named("type") = type,
		Rcpp::Named("usable") = usable,
		Rcpp::Named("stringsAsFactors") = false);
}

Rcpp::NumericMatrix distance_matrix(Rcpp::NumericMatrix X,
		const std::string& backend, int device) {
#ifdef GED_USE_WGPU
	const std::string b = normalize_backend(backend);
	if (b == "gpu" || (b == "auto" && device >= 0)) {
		return ged_wgpu_backend::distance_matrix(X, device);
	}
#endif
	return distance_matrix_cpu(X);
}

Rcpp::NumericMatrix kernel_matrix(Rcpp::NumericMatrix X,
		const std::string& kernel, int poly_s, double gamma, int device) {
#ifdef GED_USE_WGPU
	if (device >= 0) {
		return ged_wgpu_backend::kernel_matrix(X, kernel, poly_s, gamma, device);
	}
#endif
	return kernel_matrix_cpu(X, kernel, poly_s, gamma);
}

Rcpp::NumericVector objective_vals(Rcpp::NumericMatrix W,
		Rcpp::NumericMatrix Kgram, int device) {
#ifdef GED_USE_WGPU
    if (device >= 0) {
        return ged_wgpu_backend::objective_vals(W, Kgram, device);
    }
#endif
	return objective_vals_cpu(W, Kgram);
}

Rcpp::NumericVector multiple_kernel_objective_vals(Rcpp::NumericMatrix W,
        Rcpp::List Kgrams, Rcpp::NumericVector weights,
        Rcpp::NumericVector initial_objs, Rcpp::NumericVector running_sums,
        Rcpp::NumericVector max_reds, double maximum_gain_scaling, int device) {
#ifdef GED_USE_WGPU
    if (device >= 0) {
        return ged_wgpu_backend::multiple_kernel_objective_vals(W, Kgrams, weights, initial_objs, running_sums, max_reds, maximum_gain_scaling, device);
    }
#endif
    return multiple_kernel_objective_vals_cpu(W, Kgrams, weights, initial_objs, running_sums, max_reds, maximum_gain_scaling);
}

Rcpp::IntegerVector full_greedy_search(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Sinv, Rcpp::IntegerVector start_indicT, int max_iters, int device) {
#ifdef GED_USE_WGPU
    if (device >= 0) {
        return ged_wgpu_backend::full_greedy_search(X, Sinv, start_indicT, max_iters, device);
    }
#endif
    return start_indicT; // CPU dummy for now
}
Rcpp::NumericMatrix randomization_metrics(Rcpp::IntegerMatrix W, int device) {
    return randomization_metrics_cpu(W);
}

}
