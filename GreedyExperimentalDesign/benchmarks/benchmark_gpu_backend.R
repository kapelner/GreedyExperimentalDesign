#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(GreedyExperimentalDesign))

bench = function(label, f, nrep = 5L) {
	times = numeric(nrep)
	for (i in seq_len(nrep)) {
		gc()
		times[i] = system.time(f())[["elapsed"]]
	}
	cat(
		label,
		"median_sec=", round(stats::median(times), 4),
		"all_sec=", paste(round(times, 4), collapse = ","),
		"\n"
	)
	invisible(times)
}

max_abs_diff = function(x, y) {
	max(abs(x - y), na.rm = TRUE)
}

max_rel_diff = function(x, y) {
	denom = pmax(abs(x), abs(y), .Machine$double.eps)
	max(abs(x - y) / denom, na.rm = TRUE)
}

args = commandArgs(trailingOnly = TRUE)
n = if (length(args) >= 1) as.integer(args[[1]]) else 2048L
p = if (length(args) >= 2) as.integer(args[[2]]) else 16L
r = if (length(args) >= 3) as.integer(args[[3]]) else 64L
nrep = if (length(args) >= 4) as.integer(args[[4]]) else 5L
tol = if (length(args) >= 5) as.numeric(args[[5]]) else 1e-7

set.seed(123)

cat("GreedyExperimentalDesign version:", as.character(utils::packageVersion("GreedyExperimentalDesign")), "\n")
cat("Problem size: n=", n, " p=", p, " r=", r, " nrep=", nrep, "\n", sep = "")
cat("Option GreedyExperimentalDesign.use_gpu:", getOption("GreedyExperimentalDesign.use_gpu"), "\n")
cat("GPU available:\n")
print(ged_gpu_available())
cat("GPU devices:\n")
print(ged_gpu_devices())

if (!isTRUE(ged_gpu_available())) {
	stop("No usable GPU backend is available to this R session.")
}

X_check = matrix(stats::rnorm(min(n, 256L) * p), nrow = min(n, 256L))
W_check = matrix(sample(c(-1, 1), min(r, 64L) * nrow(X_check), replace = TRUE), nrow = min(r, 64L))

D_cpu = GreedyExperimentalDesign:::compute_distance_matrix_cpp(X_check)
D_gpu = compute_distance_matrix_gpu(X_check, backend = "gpu")
d_diff = max_abs_diff(D_cpu, D_gpu)
d_rel = max_rel_diff(D_cpu, D_gpu)
cat("distance max_abs_diff:", d_diff, "max_rel_diff:", d_rel, "\n")
stopifnot(d_diff <= tol || d_rel <= tol)

kernel_ids = c(poly = 1L, exponential = 2L, laplacian = 3L, inv_mult_quad = 4L, gaussian = 5L)
for (kernel_name in names(kernel_ids)) {
	K_cpu = GreedyExperimentalDesign:::compute_kernel_matrix_cpp(X_check, kernel_ids[[kernel_name]], 2L)
	K_gpu = if (kernel_name == "poly") {
		compute_kernel_matrix_gpu(X_check, kernel_name, poly_s = 2L)
	} else {
		compute_kernel_matrix_gpu(X_check, kernel_name)
	}
	k_diff = max_abs_diff(K_cpu, K_gpu)
	k_rel = max_rel_diff(K_cpu, K_gpu)
	cat("kernel", kernel_name, "max_abs_diff:", k_diff, "max_rel_diff:", k_rel, "\n")
	stopifnot(k_diff <= tol || k_rel <= tol)
}

K_check = GreedyExperimentalDesign:::compute_kernel_matrix_cpp(X_check, 5L, 2L)
O_cpu = GreedyExperimentalDesign:::compute_kernel_quadratic_forms_cpp(W_check, K_check)
O_gpu = compute_objective_vals_gpu(W_check, K_check)
o_diff = max_abs_diff(O_cpu, O_gpu)
o_rel = max_rel_diff(O_cpu, O_gpu)
cat("objective max_abs_diff:", o_diff, "max_rel_diff:", o_rel, "\n")
stopifnot(o_diff <= tol || o_rel <= tol)

X = matrix(stats::rnorm(n * p), nrow = n)
W = matrix(sample(c(-1, 1), r * n, replace = TRUE), nrow = r)
K = GreedyExperimentalDesign:::compute_kernel_matrix_cpp(X, 5L, 2L)

cat("\nBenchmarks:\n")
bench("distance_cpu", function() GreedyExperimentalDesign:::compute_distance_matrix_cpp(X), nrep)
bench("distance_gpu", function() compute_distance_matrix_gpu(X, backend = "gpu"), nrep)
bench("kernel_cpu", function() GreedyExperimentalDesign:::compute_kernel_matrix_cpp(X, 5L, 2L), nrep)
bench("kernel_gpu", function() compute_kernel_matrix_gpu(X, "gaussian"), nrep)
bench("objective_cpu", function() GreedyExperimentalDesign:::compute_kernel_quadratic_forms_cpp(W, K), nrep)
bench("objective_gpu", function() compute_objective_vals_gpu(W, K), nrep)
