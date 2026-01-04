options(java.parameters = "-Xmx10g")

SEED = 1984
N = 100
P = 5
MAX_DESIGNS = 2000
NUM_CORES = 1
PLOT_PATH = "mahal_dist_all_searches.png"
MULTI_KERNEL_NAMES = c(
  "mahalanobis",
  "gaussian",
  "laplacian",
  "exponential",
  "inv_mult_quad",
  "poly_2",
  "poly_3"
)

#!/usr/bin/env Rscript
script_path = function() {
  args = commandArgs(trailingOnly = FALSE)
  file_flag = "--file="
  match_idx = grep(file_flag, args)
  if (length(match_idx) == 0) {
    return(NULL)
  }
  normalizePath(sub(file_flag, "", args[[match_idx[1]]]))
}

this_script = script_path()
pkg_root = NULL
if (!is.null(this_script)) {
  pkg_root = normalizePath(file.path(dirname(this_script), ".."), winslash = "/", mustWork = FALSE)
}

if (!is.null(pkg_root) && requireNamespace("pkgload", quietly = TRUE) &&
  file.exists(file.path(pkg_root, "DESCRIPTION"))) {
  pkgload::load_all(pkg_root, export_all = FALSE, quiet = TRUE)
} else {
  suppressPackageStartupMessages(library(GreedyExperimentalDesign))
}

add_java_classpath = function(root_dir = NULL) {
  if (!requireNamespace("rJava", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  jar_path = NULL
  if (!is.null(root_dir)) {
    jar_path = file.path(root_dir, "inst", "java", "GreedyExperimentalDesign.jar")
    if (!file.exists(jar_path)) {
      jar_path = NULL
    }
  }
  if (is.null(jar_path)) {
    jar_path = system.file("java", "GreedyExperimentalDesign.jar", package = "GreedyExperimentalDesign")
  }
  if (!is.null(jar_path) && nzchar(jar_path) && file.exists(jar_path)) {
    rJava::.jaddClassPath(jar_path)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

add_java_classpath(pkg_root)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(MASS))
gurobi_loaded = FALSE
if (requireNamespace("gurobi", quietly = TRUE)) {
  gurobi_loaded = tryCatch({
    suppressPackageStartupMessages(library(gurobi))
    TRUE
  }, error = function(e) FALSE)
}

compute_objective_vals = getFromNamespace("compute_objective_vals_cpp", "GreedyExperimentalDesign")

set.seed(SEED)
X = MASS::mvrnorm(N, mu = rep(0, P), Sigma = diag(P))
SinvX = solve(stats::var(X))

as_integer_matrix = function(x) {
  x = as.matrix(x)
  storage.mode(x) = "integer"
  x
}

normalize_indicTs = function(indicTs, n) {
  indicTs = as.matrix(indicTs)
  if (ncol(indicTs) == n) {
    return(indicTs)
  }
  if (nrow(indicTs) == n) {
    return(t(indicTs))
  }
  stop("Unexpected indicTs dimensions: ", nrow(indicTs), "x", ncol(indicTs))
}

compute_mahal = function(indicTs, X, n, inv_cov_X) {
  indicTs = normalize_indicTs(indicTs, n)
  indicTs = as_integer_matrix(indicTs)
  compute_objective_vals(X, indicTs, "mahal_dist", inv_cov_X)
}

collect_method = function(name, indicTs, X, n, inv_cov_X) {
  vals = compute_mahal(indicTs, X, n, inv_cov_X)
  data.frame(
    method = name,
    rank = seq_along(vals),
    mahal_dist = sort(vals),
    stringsAsFactors = FALSE
  )
}

results = list()

cat("Running greedy search.\n")
ged = initGreedyExperimentalDesignObject(
  X,
  max_designs = MAX_DESIGNS,
  num_cores = NUM_CORES,
  objective = "mahal_dist",
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  verbose = FALSE
)
ged_res = resultsGreedySearch(ged, max_vectors = MAX_DESIGNS, form = "one_zero")
results[["greedy"]] = collect_method("greedy", ged_res$ending_indicTs, X, N, SinvX)

cat("Running multiple-kernel greedy search.\n")
mk = initGreedyMultipleKernelExperimentalDesignObject(
  X,
  max_designs = MAX_DESIGNS,
  kernel_pre_num_designs = MAX_DESIGNS,
  num_cores = NUM_CORES,
  kernel_names = MULTI_KERNEL_NAMES,
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  diagnostics = FALSE,
  verbose = FALSE
)
mk_res = resultsMultipleKernelGreedySearch(mk, max_vectors = MAX_DESIGNS, form = "one_zero")
results[["multiple_kernel"]] = collect_method("multiple_kernel", mk_res$ending_indicTs, X, N, SinvX)

cat("Running binary match search.\n")
bms = computeBinaryMatchStructure(X)
bm = initBinaryMatchExperimentalDesignSearchObject(
  bms,
  max_designs = MAX_DESIGNS,
  num_cores = NUM_CORES,
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  verbose = FALSE
)
bm_res = resultsBinaryMatchSearch(bm, form = "one_zero")
results[["binary_match"]] = collect_method("binary_match", bm_res, X, N, SinvX)

cat("Running binary match + greedy search.\n")
btg = initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
  X,
  max_designs = MAX_DESIGNS,
  num_cores = NUM_CORES,
  objective = "mahal_dist",
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  verbose = FALSE
)
btg_res = resultsBinaryMatchThenGreedySearch(btg, num_vectors = MAX_DESIGNS, form = "one_zero")
results[["binary_match_then_greedy"]] = collect_method("binary_match_then_greedy", btg_res$indicTs, X, N, SinvX)

cat("Running binary match + rerandomization search.\n")
btr = initBinaryMatchFollowedByRerandomizationDesignSearchObject(
  X,
  max_designs = MAX_DESIGNS * 10,
  num_cores = NUM_CORES,
  objective = "mahal_dist",
  obj_val_cutoff_to_include = Inf,
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  verbose = FALSE
)
btr_res = resultsBinaryMatchThenRerandomizationSearch(
  btr,
  num_vectors = MAX_DESIGNS * 10,
  compute_obj_vals = TRUE,
  form = "one_zero",
  use_safe_inverse = FALSE
)
btr_keep = min(length(btr_res$obj_vals), nrow(btr_res$indicTs))
btr_n_keep = max(1L, floor(btr_keep * 0.1))
btr_order = order(btr_res$obj_vals[seq_len(btr_keep)])
btr_indicTs = btr_res$indicTs[btr_order[seq_len(btr_n_keep)], , drop = FALSE]
results[["binary_match_then_rerandomization"]] = collect_method(
  "binary_match_then_rerandomization",
  btr_indicTs,
  X,
  N,
  SinvX
)

cat("Running rerandomization search.\n")
rr = initRerandomizationExperimentalDesignObject(
  X,
  max_designs = MAX_DESIGNS * 10,
  num_cores = NUM_CORES,
  objective = "mahal_dist",
  obj_val_cutoff_to_include = Inf,
  start = TRUE,
  wait = TRUE,
  seed = SEED,
  verbose = FALSE
)
rr_res = resultsRerandomizationSearch(rr, include_assignments = TRUE, form = "one_zero")
rr_vals = rr_res$obj_vals
rr_indicTs = rr_res$ending_indicTs
rr_keep = min(length(rr_vals), nrow(rr_indicTs))
rr_n_keep = max(1L, floor(rr_keep * 0.1))
rr_order = order(rr_vals[seq_len(rr_keep)])
rr_indicTs = rr_indicTs[rr_order[seq_len(rr_n_keep)], , drop = FALSE]
results[["rerandomization"]] = collect_method("rerandomization", rr_indicTs, X, N, SinvX)

cat("Running complete randomization.\n")
cr_indicTs = complete_randomization_with_forced_balanced(
  n = N,
  r = MAX_DESIGNS,
  form = "one_zero",
  seed = SEED
)
results[["complete_randomization"]] = collect_method("complete_randomization", cr_indicTs, X, N, SinvX)

if (gurobi_loaded) {
  cat("Running Gurobi search.\n")
  gurobi_results = tryCatch(
    {
      gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(
        X,
        r = MAX_DESIGNS,
        num_cores = NUM_CORES,
        objective = "mahal_dist",
        initial_time_limit_sec = 60,
        restart_time_limit_sec = 10,
        verbose = FALSE,
        use_safe_inverse = TRUE
      )
      gnoed_res = resultsGurobiNumericalOptimizeSearch(gnoed)
      collect_method("gurobi_multiple_designs", gnoed_res$indicTs, X, N, SinvX)
    },
    error = function(err) {
      cat("Skipping gurobi_multiple_designs: ", conditionMessage(err), "\n", sep = "")
      NULL
    }
  )
  if (!is.null(gurobi_results)) {
    results[["gurobi_multiple_designs"]] = gurobi_results
  }
} else {
  cat("Skipping Gurobi search: gurobi package not available.\n")
}

cat("Skipping Karp search: only supports abs_sum_diff and p = 1.\n")

cat("Preparing plot data.\n")

plot_data = do.call(rbind, results)
plot_data$method = factor(plot_data$method, levels = names(results))

method_levels = levels(plot_data$method)
mean_data = do.call(
  rbind,
  lapply(method_levels, function(method_name) {
    vals = plot_data$mahal_dist[plot_data$method == method_name]
    n = length(vals)
    sd_val = stats::sd(vals)
    se = sd_val / sqrt(n)
    ci_width = 2 * stats::qnorm(0.975) * se
    data.frame(
      method = method_name,
      mean = mean(vals),
      ci_width = max(ci_width, .Machine$double.eps),
      stringsAsFactors = FALSE
    )
  })
)

cat("Rendering plot to: ", PLOT_PATH, "\n", sep = "")

plot_obj = ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = mahal_dist, fill = method, color = method)
) +
  ggplot2::geom_histogram(alpha = 0.2, position = "identity", bins = 80, color = NA) +
  ggplot2::geom_vline(
    data = mean_data,
    ggplot2::aes(xintercept = mean, color = method, linewidth = ci_width),
    alpha = 0.6
  ) +
  ggplot2::scale_linewidth(range = c(0.6, 2.5), guide = "none") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::scale_fill_brewer(palette = "Dark2") +
  ggplot2::guides(
    fill = "none",
    color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, alpha = 0.9))
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = "Mahalanobis Distance Across Search Methods",
    x = "log10(scaled mahalanobis distance)",
    y = "",
    color = "Method"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(PLOT_PATH, plot_obj, width = 9, height = 6, dpi = 120)
cat("Saved plot to:", PLOT_PATH, "\n")
invisible(NULL)
