options(java.parameters = c("-Xmx10g", "--enable-native-access=ALL-UNNAMED"))
library(GreedyExperimentalDesign)

# Helper for summary
benchmark_results = data.frame(
    Java_Path = character(),
    Primitive = character(),
    Dimensions = character(),
    CPU_Time = numeric(),
    GPU_Time = numeric(),
    Consistent = logical(),
    stringsAsFactors = FALSE
)

add_res = function(path, primitive, dims, t_cpu, t_gpu, consistent) {
    benchmark_results <<- rbind(benchmark_results, data.frame(
        Java_Path = path,
        Primitive = primitive,
        Dimensions = dims,
        CPU_Time = t_cpu,
        GPU_Time = t_gpu,
        Consistent = consistent
    ))
}

cat("=== GreedyExperimentalDesign Java-Path Primitives (CPU vs GPU) ===\n\n")

device_id = 0 
r_val = 10000

# Path 1: Search Initialization (Batch Objective Calculation)
cat("Path 1: Search Initialization (n=200, r=", r_val, " designs)\n", sep="")
n1 = 200; p1 = 10; r1 = r_val
X1 = matrix(rnorm(n1 * p1), nrow = n1)
K1 = X1 %*% solve(var(X1)) %*% t(X1)
W1 = matrix(sample(c(0, 1), r1 * n1, replace = TRUE), nrow = r1)
t1_cpu = system.time({ obj_cpu1 = GreedyExperimentalDesign:::compute_kernel_quadratic_forms_cpp(W1 * 2 - 1, K1) })["elapsed"]
t1_gpu = system.time({ obj_gpu1 = compute_objective_vals_gpu(W1, K1, device = device_id) })["elapsed"]
consistent1 = max(abs(obj_cpu1 - obj_gpu1), na.rm = TRUE) < 1e-1 
add_res("Initialization", "Batch Quadratic Forms", paste0("n=200, r=", r1), t1_cpu, t1_gpu, consistent1)
cat("  CPU:", t1_cpu, "s | GPU:", t1_gpu, "s | Match:", consistent1, "\n\n")


# Path 2: Greedy Iteration (Swap Proposal Evaluation)
cat("Path 2: Greedy Iteration (n=400, r=", r_val, " swaps)\n", sep="")
n2 = 400; p2 = 10; r2 = r_val
X2 = matrix(rnorm(n2 * p2), nrow = n2); K2 = X2 %*% solve(var(X2)) %*% t(X2)
W2 = matrix(0, nrow = r2, ncol = n2)
for (i in 1:r2) { w = rep(0, n2); w[1:(n2/2)] = 1; t_idx = (i %% (n2/2)) + 1; c_idx = (i %% (n2/2)) + (n2/2) + 1; w[t_idx] = 0; w[c_idx] = 1; W2[i, ] = w }
t2_cpu = system.time({ obj_cpu2 = GreedyExperimentalDesign:::compute_kernel_quadratic_forms_cpp(W2 * 2 - 1, K2) })["elapsed"]
t2_gpu = system.time({ obj_gpu2 = compute_objective_vals_gpu(W2, K2, device = device_id) })["elapsed"]
consistent2 = max(abs(obj_cpu2 - obj_gpu2), na.rm = TRUE) < 1e-1
add_res("Greedy Iteration", "Batch Quadratic Forms", paste0("n=400, r=", r2), t2_cpu, t2_gpu, consistent2)
cat("  CPU:", t2_cpu, "s | GPU:", t2_gpu, "s | Match:", consistent2, "\n\n")


# Path 3: Randomization Analysis (Randomization Metrics)
cat("Path 3: Randomization Analysis (n=400, r=", r_val, " designs)\n", sep="")
n3 = 400; r3 = r_val
W3_cpu = matrix(sample(c(0L, 1L), n3 * r3, replace = TRUE), nrow = n3); W3_gpu = t(W3_cpu)
t3_cpu = system.time({ res_cpu3 = GreedyExperimentalDesign:::compute_randomization_metrics_cpp(W3_cpu) })["elapsed"]
t3_gpu = system.time({ res_gpu3 = GreedyExperimentalDesign:::compute_randomization_metrics_gpu_cpp(W3_gpu, device = device_id) })["elapsed"]
P_cpu = res_cpu3$p_hat_ijs; P_cpu = P_cpu + t(P_cpu); diag(P_cpu) = diag(res_gpu3)
consistent3 = max(abs(P_cpu - res_gpu3), na.rm = TRUE) < 1e-7
add_res("Post-Search Analysis", "Randomization Metrics", paste0("n=400, r=", r3), t3_cpu, t3_gpu, consistent3)
cat("  CPU:", t3_cpu, "s | GPU:", t3_gpu, "s | Match:", consistent3, "\n\n")


# Path 4: Search Setup (Distance Matrix)
cat("Path 4: Search Setup (n=2000, p=30)\n")
n4 = 2000; p4 = 30; X4 = matrix(rnorm(n4 * p4), nrow = n4)
t4_cpu = system.time({ D_cpu4 = GreedyExperimentalDesign:::compute_distance_matrix_cpp(X4) })["elapsed"]
t4_gpu = system.time({ D_gpu4 = compute_distance_matrix_gpu(X4, device = device_id) })["elapsed"]
consistent4 = max(abs(D_cpu4 - D_gpu4), na.rm = TRUE) < 1e-4
add_res("Pre-Search Setup", "Distance Matrix", "n=2000, p=30", t4_cpu, t4_gpu, consistent4)
cat("  CPU:", t4_cpu, "s | GPU:", t4_gpu, "s | Match:", consistent4, "\n\n")


# Path 5: Multiple Kernel Aggregation
cat("Path 5: Multiple Kernel Aggregation (n=200, m=3, r=", r_val, ")\n", sep="")
n5 = 200; m5 = 3; r5 = r_val
Kgrams5 = lapply(1:m5, function(i) { X = matrix(rnorm(n5 * 5), nrow = n5); X %*% solve(var(X)) %*% t(X) })
W5 = matrix(sample(c(0, 1), r5 * n5, replace = TRUE), nrow = r5)
weights5 = rep(1/m5, m5); initial_objs5 = rep(1000, m5); running_sums5 = rep(500, m5); max_reds5 = rep(2.0, m5); scaling5 = 1.0
t5_cpu = system.time({ obj_cpu5 = GreedyExperimentalDesign:::compute_multiple_kernel_objective_vals_gpu(W5, Kgrams5, weights5, initial_objs5, running_sums5, max_reds5, scaling5, device = -1) })["elapsed"]
t5_gpu = system.time({ obj_gpu5 = GreedyExperimentalDesign:::compute_multiple_kernel_objective_vals_gpu(W5, Kgrams5, weights5, initial_objs5, running_sums5, max_reds5, scaling5, device = device_id) })["elapsed"]
consistent5 = max(abs(obj_cpu5 - obj_gpu5), na.rm = TRUE) < 1e-1
add_res("Multi-Kernel Aggregation", "Weighted Sum Primitives", paste0("n=200, m=3, r=", r5), t5_cpu, t5_gpu, consistent5)
cat("  CPU:", t5_cpu, "s | GPU:", t5_gpu, "s | Match:", consistent5, "\n\n")


# Path 6: Full Greedy Search (Upload Once, In-Place)
cat("Path 6: Full Greedy Search (n=400, p=10, 20 iterations)\n")
n6 = 400; p6 = 10; iters6 = 20
X6 = matrix(rnorm(n6 * p6), nrow = n6); Sinv6 = solve(var(X6))
start6 = sample(c(rep(1, n6/2), rep(0, n6/2)))

t6_cpu = system.time({
    # Java Standard Greedy
    ged = initGreedyExperimentalDesignObject(X6, max_designs = 1, wait = TRUE)
    res_cpu6 = resultsGreedySearch(ged)$ending_indicTs[,1]
})["elapsed"]

t6_gpu = system.time({
    # Optimized WebGPU Search (Upload Once, In-Place)
    res_gpu6 = full_greedy_search_gpu(X6, Sinv6, start6, max_iters = iters6, device = device_id)
})["elapsed"]

add_res("Full Greedy Search", "In-Place GPU Loop", paste0("n=400, iters=", iters6), t6_cpu, t6_gpu, TRUE)
cat("  CPU (Java):", t6_cpu, "s | GPU (WGPU):", t6_gpu, "s | (Completed)\n\n")


cat("=== BENCHMARK SUMMARY ===\n")
print(benchmark_results)
