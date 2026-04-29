## Run with: Rscript GreedyExperimentalDesign/benchmarks/run_version_benchmarks.R --lib <lib> --output <csv> --label <label>
JAVA_XMS_GB = 10
JAVA_XMX_GB = 10
MIN_SECONDS = 1
MAX_INNER_REPEATS = 200
MIN_ELAPSED_SEC = 1e-06
TIME_REPS = 3
TIME_WARMUP = 1
SEED = 1984
NUM_CORES = 1
MAX_DESIGNS = 10
N = 40
P = 2
FIRST_COL = 1
DIAGNOSTICS = FALSE
DIAGNOSTICS_MULTIPLE_KERNEL = FALSE
VERBOSE = FALSE
DIFF_METHOD = FALSE
USE_SAFE_INVERSE = FALSE

get_arg = function(flag, default = NULL) {
  args = commandArgs(trailingOnly = TRUE)
  idx = match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

has_flag = function(flag) {
  args = commandArgs(trailingOnly = TRUE)
  flag %in% args
}

require_tool = function(tool) {
  tool_path = Sys.which(tool)
  if (tool_path == "") {
    stop("Required tool not found on PATH: ", tool)
  }
  tool_path
}

run_ant_dist = function(repo_dir, ant_bin) {
  build_file = file.path(repo_dir, "build.xml")
  if (!file.exists(build_file)) {
    stop("build.xml not found in repo: ", repo_dir)
  }
  old_wd = getwd()
  setwd(repo_dir)
  on.exit(setwd(old_wd), add = TRUE)
  status = system2(ant_bin, "dist")
  if (status != 0) {
    stop("Ant build failed in repo: ", repo_dir)
  }
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

lib_path = get_arg("--lib")
out_path = get_arg("--output")
repo_arg = get_arg("--repo")
label = get_arg("--label", "run")
quiet = has_flag("--quiet")
if (is.null(lib_path) || lib_path == "") {
  stop("Missing --lib argument.")
}
if (is.null(out_path) || out_path == "") {
  stop("Missing --output argument.")
}
if (!dir.exists(lib_path)) {
  stop("Library path does not exist: ", lib_path)
}

options(java.parameters = sprintf("-Xmx%sg", JAVA_XMX_GB))

script_path = function() {
  args = commandArgs(trailingOnly = FALSE)
  file_flag = "--file="
  match_idx = grep(file_flag, args)
  if (length(match_idx) == 0) {
    return(NULL)
  }
  normalizePath(sub(file_flag, "", args[[match_idx[1]]]))
}

ant_bin = require_tool("ant")
if (!is.null(repo_arg)) {
  pkg_dir = normalizePath(repo_arg)
} else {
  this_script = script_path()
  if (is.null(this_script)) {
    pkg_dir = normalizePath(getwd())
  } else {
    pkg_dir = normalizePath(file.path(dirname(this_script), "..", ".."))
  }
}
cat("Rebuilding Java jar in:", pkg_dir, "\n")
run_ant_dist(pkg_dir, ant_bin)

suppressPackageStartupMessages(library(GreedyExperimentalDesign, lib.loc = lib_path))

call_init = function(fun, args, verbose = VERBOSE, suppress = FALSE) {
  formals_names = names(formals(fun))
  if (!is.null(formals_names) && "verbose" %in% formals_names) {
    args$verbose = verbose
  } else if (!is.null(formals_names) && "compute_dist_matrix" %in% formals_names) {
    fun_name = deparse(substitute(fun))
    needs_compat = fun_name %in% c(
      "initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject",
      "initBinaryMatchFollowedByRerandomizationDesignSearchObject"
    )
    if (needs_compat && is.null(args$compute_dist_matrix)) {
      args$compute_dist_matrix = FALSE
    }
  }
  if (!is.null(formals_names) && !"..." %in% formals_names) {
    arg_names = names(args)
    if (!is.null(arg_names)) {
      keep = (arg_names == "" | is.na(arg_names) | arg_names %in% formals_names)
      args = args[keep]
    }
  }
  if (suppress) {
    suppressWarnings(do.call(fun, args))
  } else {
    do.call(fun, args)
  }
}

measure_elapsed = function(fun) {
  start = proc.time()[["elapsed"]]
  fun()
  elapsed = proc.time()[["elapsed"]] - start
  if (!is.finite(elapsed) || elapsed <= 0) {
    elapsed = MIN_ELAPSED_SEC
  }
  elapsed
}

calc_inner_repeats = function(elapsed, min_seconds = MIN_SECONDS, max_repeats = MAX_INNER_REPEATS) {
  repeats = ceiling(min_seconds / elapsed)
  repeats = max(1, repeats)
  min(repeats, max_repeats)
}

time_expr = function(fun, reps = TIME_REPS, warmup = TIME_WARMUP, min_seconds = MIN_SECONDS,
                     max_repeats = MAX_INNER_REPEATS) {
  stopifnot(reps >= 1, warmup >= 0)
  gc()
  for (i in seq_len(warmup)) {
    fun()
  }
  gc()
  inner_repeats = calc_inner_repeats(measure_elapsed(fun), min_seconds, max_repeats)
  times = numeric(reps)
  for (i in seq_len(reps)) {
    start = proc.time()[["elapsed"]]
    for (j in seq_len(inner_repeats)) {
      fun()
    }
    times[i] = (proc.time()[["elapsed"]] - start) / inner_repeats
  }
  list(times = times, inner_repeats = inner_repeats)
}

time_expr_with_setup = function(setup_fun, run_fun, reps = TIME_REPS, warmup = TIME_WARMUP,
                                min_seconds = MIN_SECONDS, max_repeats = MAX_INNER_REPEATS) {
  stopifnot(reps >= 1, warmup >= 0)
  gc()
  for (i in seq_len(warmup)) {
    obj = setup_fun()
    run_fun(obj)
  }
  gc()
  obj = setup_fun()
  inner_repeats = calc_inner_repeats(measure_elapsed(function() run_fun(obj)), min_seconds, max_repeats)
  rm(obj)
  gc()
  times = numeric(reps)
  for (i in seq_len(reps)) {
    if (inner_repeats == 1) {
      obj = setup_fun()
      start = proc.time()[["elapsed"]]
      run_fun(obj)
      times[i] = proc.time()[["elapsed"]] - start
      rm(obj)
      gc()
    } else {
      objs = vector("list", inner_repeats)
      for (j in seq_len(inner_repeats)) {
        objs[[j]] = setup_fun()
      }
      start = proc.time()[["elapsed"]]
      for (j in seq_len(inner_repeats)) {
        run_fun(objs[[j]])
      }
      times[i] = (proc.time()[["elapsed"]] - start) / inner_repeats
      rm(objs)
      gc()
    }
  }
  list(times = times, inner_repeats = inner_repeats)
}

bench_case = function(name, fun, reps = TIME_REPS, warmup = TIME_WARMUP) {
  timing = time_expr(fun, reps = reps, warmup = warmup)
  data.frame(
    name = name,
    rep = seq_along(timing$times),
    inner_repeats = rep(timing$inner_repeats, length(timing$times)),
    elapsed_sec = timing$times,
    stringsAsFactors = FALSE
  )
}

bench_case_with_setup = function(name, setup_fun, run_fun, reps = TIME_REPS, warmup = TIME_WARMUP) {
  timing = time_expr_with_setup(setup_fun, run_fun, reps = reps, warmup = warmup)
  data.frame(
    name = name,
    rep = seq_along(timing$times),
    inner_repeats = rep(timing$inner_repeats, length(timing$times)),
    elapsed_sec = timing$times,
    stringsAsFactors = FALSE
  )
}

set.seed(SEED)
X = matrix(rnorm(N * P), nrow = N, ncol = P)
X1 = X[, FIRST_COL, drop = FALSE]

make_greedy = function(start, wait) {
  call_init(
    initGreedyExperimentalDesignObject,
    list(
      X,
      max_designs = MAX_DESIGNS,
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED,
      diagnostics = DIAGNOSTICS
    )
  )
}

make_multiple_kernel = function(start, wait) {
  call_init(
    initGreedyMultipleKernelExperimentalDesignObject,
    list(
      X,
      max_designs = MAX_DESIGNS,
      kernel_pre_num_designs = MAX_DESIGNS,
      kernel_names = c("mahalanobis", "gaussian"),
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED,
      diagnostics = DIAGNOSTICS_MULTIPLE_KERNEL,
      use_safe_inverse = USE_SAFE_INVERSE
    ),
    suppress = TRUE
  )
}

make_binary_match = function(start, wait) {
  bms = computeBinaryMatchStructure(X)
  call_init(
    initBinaryMatchExperimentalDesignSearchObject,
    list(
      bms,
      max_designs = MAX_DESIGNS,
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED
    )
  )
}

make_binary_match_then_greedy = function(start, wait) {
  call_init(
    initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject,
    list(
      X,
      diff_method = DIFF_METHOD,
      max_designs = MAX_DESIGNS,
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED,
      diagnostics = DIAGNOSTICS
    )
  )
}

make_binary_match_then_rerand = function(start, wait) {
  call_init(
    initBinaryMatchFollowedByRerandomizationDesignSearchObject,
    list(
      X,
      max_designs = MAX_DESIGNS,
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED,
      obj_val_cutoff_to_include = Inf
    )
  )
}

make_rerandomization = function(start, wait) {
  call_init(
    initRerandomizationExperimentalDesignObject,
    list(
      X,
      obj_val_cutoff_to_include = Inf,
      max_designs = MAX_DESIGNS,
      wait = wait,
      start = start,
      num_cores = NUM_CORES,
      seed = SEED
    )
  )
}

make_karp = function(start, wait) {
  call_init(
    initKarpExperimentalDesignObject,
    list(
      X1,
      wait = wait,
      balanced = TRUE,
      start = start
    )
  )
}

bench_specs = list(
  list(
    name = "greedy_init_only",
    fun = function() make_greedy(start = FALSE, wait = FALSE)
  ),
  list(
    name = "greedy_search_only",
    setup = function() make_greedy(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj)
  ),
  list(
    name = "greedy_full_search",
    fun = function() make_greedy(start = TRUE, wait = TRUE)
  ),
  list(
    name = "multiple_kernel_init_only",
    fun = function() make_multiple_kernel(start = FALSE, wait = FALSE)
  ),
  list(
    name = "multiple_kernel_search_only",
    setup = function() make_multiple_kernel(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj)
  ),
  list(
    name = "multiple_kernel_full_search",
    fun = function() make_multiple_kernel(start = TRUE, wait = TRUE)
  ),
  list(
    name = "binary_match_init_only",
    fun = function() make_binary_match(start = FALSE, wait = FALSE)
  ),
  list(
    name = "binary_match_search_only",
    setup = function() make_binary_match(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj)
  ),
  list(
    name = "binary_match_full_search",
    fun = function() make_binary_match(start = TRUE, wait = TRUE)
  ),
  list(
    name = "binary_match_then_greedy_init_only",
    fun = function() make_binary_match_then_greedy(start = FALSE, wait = FALSE)
  ),
  list(
    name = "binary_match_then_greedy_search_only",
    setup = function() make_binary_match_then_greedy(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj$greedy_design)
  ),
  list(
    name = "binary_match_then_greedy_full_search",
    fun = function() make_binary_match_then_greedy(start = TRUE, wait = TRUE)
  ),
  list(
    name = "binary_match_then_rerandomization_init_only",
    fun = function() make_binary_match_then_rerand(start = FALSE, wait = FALSE)
  ),
  list(
    name = "binary_match_then_rerandomization_search_only",
    setup = function() make_binary_match_then_rerand(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj$rerandomization_design)
  ),
  list(
    name = "binary_match_then_rerandomization_full_search",
    fun = function() make_binary_match_then_rerand(start = TRUE, wait = TRUE)
  ),
  list(
    name = "rerandomization_init_only",
    fun = function() make_rerandomization(start = FALSE, wait = FALSE)
  ),
  list(
    name = "rerandomization_search_only",
    setup = function() make_rerandomization(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj)
  ),
  list(
    name = "rerandomization_full_search",
    fun = function() make_rerandomization(start = TRUE, wait = TRUE)
  ),
  list(
    name = "karp_init_only",
    fun = function() make_karp(start = FALSE, wait = FALSE)
  ),
  list(
    name = "karp_search_only",
    setup = function() make_karp(start = FALSE, wait = TRUE),
    run = function(obj) startSearch(obj)
  ),
  list(
    name = "karp_full_search",
    fun = function() make_karp(start = TRUE, wait = TRUE)
  )
)

results_list = vector("list", length(bench_specs))
for (i in seq_along(bench_specs)) {
  spec = bench_specs[[i]]
  if (!quiet) {
    cat("Running benchmark:", spec$name, "\n")
  }
  if (!is.null(spec$setup)) {
    results_list[[i]] = bench_case_with_setup(spec$name, spec$setup, spec$run)
  } else {
    results_list[[i]] = bench_case(spec$name, spec$fun)
  }
}

results = do.call(rbind, results_list)
results$label = label
results$version = as.character(utils::packageVersion("GreedyExperimentalDesign"))

write.csv(results, out_path, row.names = FALSE)
if (!quiet) {
  print(results, row.names = FALSE)
}
