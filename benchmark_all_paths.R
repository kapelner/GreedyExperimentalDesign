library(GreedyExperimentalDesign)

# Helper for summary
benchmark_results = data.frame(
    Path = character(),
    Size = character(),
    Native_CPP_Time = numeric(),
    WebGPU_Time = numeric(),
    Speedup = character(),
    Consistent = logical(),
    stringsAsFactors = FALSE
)

add_res = function(path, size, t_cpp, t_gpu, consistent) {
    speedup = paste0(round(t_cpp / t_gpu, 2), "x")
    benchmark_results <<- rbind(benchmark_results, data.frame(
        Path = path,
        Size = size,
        Native_CPP_Time = t_cpp,
        WebGPU_Time = t_gpu,
        Speedup = speedup,
        Consistent = consistent
    ))
}

cat("=== GreedyExperimentalDesign FULL BACKEND BENCHMARK (C++ vs WebGPU) ===\n\n")

# Path 1: Exhaustive Search (n=22)
cat("Path 1: Exhaustive Search (n=22, p=5, r=705432 designs)\n")
n1 = 22; p1 = 5
X1 = matrix(rnorm(n1 * p1), nrow = n1)
Sinv1 = solve(var(X1))
K1 = X1 %*% Sinv1 %*% t(X1)
combos1 = combn(n1, n1/2)
W1 = matrix(0, nrow = ncol(combos1), ncol = n1)
for (i in 1:ncol(combos1)) { W1[i, combos1[,i]] = 1 }

t1_cpp = system.time({
    obj_cpp1 = apply(W1, 1, function(w) {
        v = 2*w - 1
        as.numeric(t(v) %*% K1 %*% v)
    })
})["elapsed"]

t1_gpu = system.time({
    obj_gpu1 = compute_objective_vals_gpu(W1, K1)
})["elapsed"]

diff1 = max(abs(obj_cpp1 - obj_gpu1))
cat("  Max Diff Path 1:", diff1, "\n")
consistent1 = diff1 < 1e-3
add_res("Exhaustive Search", "n=22, r=705k", t1_cpp, t1_gpu, consistent1)
cat("  C++:", t1_cpp, "s | WebGPU:", t1_gpu, "s | Match:", consistent1, "\n\n")


# Path 2: Greedy Search Proposals (n=200)
cat("Path 2: Greedy Search Proposals (n=200, p=10, 10000 swaps)\n")
n2 = 200; p2 = 10
X2 = matrix(rnorm(n2 * p2), nrow = n2)
Sinv2 = solve(var(X2))
K2 = X2 %*% Sinv2 %*% t(X2)
W2 = matrix(0, nrow = 10000, ncol = n2)
for(i in 1:10000) W2[i, sample(1:n2, n2/2)] = 1

t2_cpp = system.time({
    obj_cpp2 = apply(W2, 1, function(w) {
        v = 2*w - 1
        as.numeric(t(v) %*% K2 %*% v)
    })
})["elapsed"]

t2_gpu = system.time({
    obj_gpu2 = compute_objective_vals_gpu(W2, K2)
})["elapsed"]

diff2 = max(abs(obj_cpp2 - obj_gpu2))
cat("  Max Diff Path 2:", diff2, "\n")
consistent2 = diff2 < 0.005
add_res("Greedy Proposals", "n=200, r=10k", t2_cpp, t2_gpu, consistent2)
cat("  C++:", t2_cpp, "s | WebGPU:", t2_gpu, "s | Match:", consistent2, "\n\n")


# Path 3: Randomization Metrics (n=32, r=1000)
cat("Path 3: Randomization Metrics (n=32, r=1000 designs)\n")
n3 = 32; r3 = 1000
W3 = matrix(sample(c(0L, 1L), r3 * n3, replace = TRUE), nrow = r3)

t3_cpp = system.time({
    P_cpp = compute_randomization_metrics_gpu(W3, device = -1) 
})["elapsed"]

t3_gpu = system.time({
    P_gpu = compute_randomization_metrics_gpu(W3)
})["elapsed"]

consistent3 = max(abs(P_cpp - P_gpu), na.rm = TRUE) < 1e-7
add_res("Randomization Metrics", "n=32, r=1000", t3_cpp, t3_gpu, consistent3)
cat("  C++:", t3_cpp, "s | WebGPU:", t3_gpu, "s | Match:", consistent3, "\n\n")


# Path 4: Kernel Matrix (n=2000)
cat("Path 4: Kernel Matrix (n=2000, p=10)\n")
n4 = 2000; p4 = 10
X4 = matrix(rnorm(n4 * p4), nrow = n4)

t4_cpp = system.time({
    K_cpp = GreedyExperimentalDesign:::compute_kernel_matrix_cpp(X4, 5, 2)
})["elapsed"]

t4_gpu = system.time({
    K_gpu = compute_kernel_matrix_gpu(X4, kernel = "gaussian")
})["elapsed"]

consistent4 = max(abs(K_cpp - K_gpu)) < 1e-7
add_res("Kernel Matrix", "n=2000", t4_cpp, t4_gpu, consistent4)
cat("  C++:", t4_cpp, "s | WebGPU:", t4_gpu, "s | Match:", consistent4, "\n\n")

cat("=== BENCHMARK SUMMARY ===\n")
print(benchmark_results)
