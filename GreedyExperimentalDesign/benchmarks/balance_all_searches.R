#!/usr/bin/env Rscript
## Balance metric comparison across all main design paths.
## Usage: Rscript balance_all_searches.R <lib_path> <label> [pregpu]
args = commandArgs(trailingOnly = TRUE)
lib_path  = args[[1]]
label     = args[[2]]
is_pregpu = length(args) >= 3 && args[[3]] == "pregpu"

options(java.parameters = "-Xmx6g")
suppressPackageStartupMessages(library(GreedyExperimentalDesign, lib.loc = lib_path))

SEED        = 1984
N           = 100
P           = 5
MAX_DESIGNS = 500
NUM_CORES   = 1
MK_KERNELS  = c("mahalanobis", "gaussian", "laplacian", "exponential")

set.seed(SEED)
X     = MASS::mvrnorm(N, mu = rep(0, P), Sigma = diag(P))
SinvX = solve(var(X))

## Helper: Mahal distance for each row of a 0/1 indicator matrix (r x n)
mahal_vals = function(indicTs) {
  if (is.null(indicTs) || length(indicTs) == 0) return(numeric(0))
  if (is.vector(indicTs)) indicTs = matrix(indicTs, nrow = 1)
  apply(indicTs, 1, function(w) {
    iT = which(w == 1); iC = which(w == 0)
    d  = colMeans(X[iT, , drop = FALSE]) - colMeans(X[iC, , drop = FALSE])
    as.numeric(t(d) %*% SinvX %*% d)
  })
}

## Compat wrappers: pre-GPU lacks 'verbose'; some args differ
init_greedy = function(...) {
  args_list = list(...)
  if (is_pregpu) args_list$verbose = NULL
  do.call(initGreedyExperimentalDesignObject, args_list)
}
init_rerand = function(...) {
  args_list = list(...)
  if (is_pregpu) args_list$verbose = NULL
  do.call(initRerandomizationExperimentalDesignObject, args_list)
}
init_bm = function(...) {
  args_list = list(...)
  if (is_pregpu) {
    args_list$verbose = NULL
    do.call(initBinaryMatchExperimentalDesignSearch, args_list)
  } else {
    do.call(initBinaryMatchExperimentalDesignSearchObject, args_list)
  }
}
init_btg = function(...) {
  args_list = list(...)
  if (is_pregpu) {
    args_list$verbose = NULL
    do.call(initBinaryMatchFollowedByGreedyExperimentalDesignSearch, args_list)
  } else {
    do.call(initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject, args_list)
  }
}
init_btr = function(...) {
  args_list = list(...)
  if (is_pregpu) {
    args_list$verbose = NULL
    do.call(initBinaryMatchFollowedByRerandomizationDesignSearch, args_list)
  } else {
    do.call(initBinaryMatchFollowedByRerandomizationDesignSearchObject, args_list)
  }
}

results = list()

## ── 1. Complete randomization (baseline) ────────────────────────────────────
cat("complete_randomization...\n")
set.seed(SEED)
half = N / 2L
cr   = t(replicate(MAX_DESIGNS, sample(c(rep(1L, half), rep(0L, half)))))
results[["complete_randomization"]] = sort(mahal_vals(cr))

## ── 2. Greedy search (mahal_dist) ───────────────────────────────────────────
cat("greedy_mahal...\n")
ged = init_greedy(X = X, objective = "mahal_dist",
                  max_designs = MAX_DESIGNS, num_cores = NUM_CORES,
                  seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE)
gr  = resultsGreedySearch(ged, max_vectors = MAX_DESIGNS, form = "one_zero")
results[["greedy_mahal"]] = sort(mahal_vals(gr$ending_indicTs))

## ── 3. Greedy search (abs_sum_diff) ─────────────────────────────────────────
cat("greedy_abs...\n")
gea = init_greedy(X = X, objective = "abs_sum_diff",
                  max_designs = MAX_DESIGNS, num_cores = NUM_CORES,
                  seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE)
gra = resultsGreedySearch(gea, max_vectors = MAX_DESIGNS, form = "one_zero")
results[["greedy_abs"]] = sort(mahal_vals(gra$ending_indicTs))

## ── 4. Rerandomization ──────────────────────────────────────────────────────
cat("rerandomization...\n")
rr = init_rerand(X = X, max_designs = MAX_DESIGNS * 10, num_cores = NUM_CORES,
                 objective = "mahal_dist", obj_val_cutoff_to_include = Inf,
                 seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE)
rr_res  = resultsRerandomizationSearch(rr, include_assignments = TRUE, form = "one_zero")
rr_k    = min(length(rr_res$obj_vals), nrow(rr_res$ending_indicTs))
rr_ord  = order(rr_res$obj_vals[seq_len(rr_k)])[seq_len(max(1L, floor(rr_k * 0.1)))]
results[["rerandomization"]] = sort(mahal_vals(rr_res$ending_indicTs[rr_ord, , drop = FALSE]))

## ── 5. Binary match ─────────────────────────────────────────────────────────
cat("binary_match...\n")
bms = computeBinaryMatchStructure(X)
bm  = init_bm(bms, max_designs = MAX_DESIGNS, num_cores = NUM_CORES,
               seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE)
bm_res = resultsBinaryMatchSearch(bm, form = "one_zero")
results[["binary_match"]] = sort(mahal_vals(bm_res))

## ── 6. Binary match → greedy ────────────────────────────────────────────────
cat("binary_match_then_greedy...\n")
btg_extra = if (is_pregpu) list(compute_dist_matrix = FALSE) else list()
btg = do.call(init_btg, c(list(X = X, max_designs = MAX_DESIGNS, num_cores = NUM_CORES,
               objective = "mahal_dist", seed = SEED,
               start = TRUE, wait = TRUE, verbose = FALSE), btg_extra))
btg_res = resultsBinaryMatchThenGreedySearch(btg, num_vectors = MAX_DESIGNS, form = "one_zero")
results[["binary_match_then_greedy"]] = sort(mahal_vals(btg_res$indicTs))

## ── 7. Binary match → rerandomization ──────────────────────────────────────
cat("binary_match_then_rerandomization...\n")
if (is_pregpu) {
  results[["binary_match_then_rerandomization"]] = NA_real_
} else {
  btr = init_btr(X = X, max_designs = MAX_DESIGNS * 10, num_cores = NUM_CORES,
                 objective = "mahal_dist", obj_val_cutoff_to_include = Inf,
                 seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE)
  btr_res = resultsBinaryMatchThenRerandomizationSearch(
    btr, num_vectors = MAX_DESIGNS * 10,
    compute_obj_vals = TRUE, form = "one_zero", use_safe_inverse = FALSE
  )
  btr_k   = min(length(btr_res$obj_vals), nrow(btr_res$indicTs))
  btr_ord = order(btr_res$obj_vals[seq_len(btr_k)])[seq_len(max(1L, floor(btr_k * 0.1)))]
  results[["binary_match_then_rerandomization"]] = sort(mahal_vals(btr_res$indicTs[btr_ord, , drop = FALSE]))
}

## ── 8. Multiple kernel greedy search ────────────────────────────────────────
cat("multiple_kernel...\n")
mk_res_indicTs = tryCatch({
  if (is_pregpu) {
    mk = suppressWarnings(initGreedyMultipleKernelExperimentalDesignObject(
      X = X, max_designs = MAX_DESIGNS, kernel_pre_num_designs = MAX_DESIGNS,
      kernel_names = MK_KERNELS, maximum_gain_scaling = 1.1,
      num_cores = NUM_CORES, seed = SEED, start = TRUE, wait = TRUE
    ))
  } else {
    mk = suppressWarnings(initGreedyMultipleKernelExperimentalDesignObject(
      X = X, max_designs = MAX_DESIGNS, kernel_pre_num_designs = MAX_DESIGNS,
      kernel_names = MK_KERNELS, num_cores = NUM_CORES,
      seed = SEED, start = TRUE, wait = TRUE, verbose = FALSE
    ))
  }
  mr = resultsMultipleKernelGreedySearch(mk, max_vectors = MAX_DESIGNS, form = "one_zero")
  mr$ending_indicTs
}, error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
results[["multiple_kernel"]] = if (!is.null(mk_res_indicTs)) sort(mahal_vals(mk_res_indicTs)) else NA_real_

## ── Print & save ────────────────────────────────────────────────────────────
rows = lapply(names(results), function(nm) {
  v = results[[nm]]
  v = v[is.finite(v)]
  if (length(v) == 0)
    return(data.frame(method = nm, label = label, n = 0L,
                      min = NA_real_, pct10 = NA_real_, median = NA_real_, mean = NA_real_))
  data.frame(method = nm, label = label, n = length(v),
             min    = min(v),
             pct10  = unname(quantile(v, 0.1)),
             median = median(v),
             mean   = mean(v))
})
out = do.call(rbind, rows)
rownames(out) = NULL

cat("\n=== Balance metrics (Mahal distance — lower is better) ===\n")
cat(sprintf("  Label: %s\n\n", label))
for (i in seq_len(nrow(out))) {
  r = out[i, ]
  cat(sprintf("  %-38s  min=%7.4f  p10=%7.4f  median=%7.4f  mean=%7.4f  (n=%d)\n",
              r$method, r$min, r$pct10, r$median, r$mean, r$n))
}
write.csv(out, paste0("/tmp/ged_balance_", label, ".csv"), row.names = FALSE)
cat(sprintf("\nSaved to /tmp/ged_balance_%s.csv\n", label))
