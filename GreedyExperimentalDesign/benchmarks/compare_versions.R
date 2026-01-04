## Run with: Rscript GreedyExperimentalDesign/benchmarks/compare_versions.R [--commit <hash>] [--output <csv>] [--keep-temp]
PREV_COMMIT = "f2433fd3555c46181c3d20a32c9aaae241d87292"

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

script_path = function() {
  args = commandArgs(trailingOnly = FALSE)
  file_flag = "--file="
  match_idx = grep(file_flag, args)
  if (length(match_idx) == 0) {
    return(NULL)
  }
  normalizePath(sub(file_flag, "", args[[match_idx[1]]]))
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

install_pkg = function(pkg_dir, lib_dir, r_bin) {
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  status = system2(r_bin, c("CMD", "INSTALL", "-l", lib_dir, pkg_dir))
  if (status != 0) {
    stop("R CMD INSTALL failed for: ", pkg_dir)
  }
  pkg_name = basename(pkg_dir)
  if (!dir.exists(file.path(lib_dir, pkg_name))) {
    stop("Package was not installed in library: ", lib_dir)
  }
}

run_bench = function(rscript_bin, script_path, lib_dir, output_path, label) {
  args = c(
    "--vanilla",
    script_path,
    "--lib", lib_dir,
    "--output", output_path,
    "--label", label,
    "--quiet"
  )
  status = system2(rscript_bin, args)
  if (status != 0) {
    stop("Benchmark run failed for label: ", label)
  }
}

summarize_runs = function(df) {
  stats = aggregate(
    elapsed_sec ~ name + label + version,
    data = df,
    FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x))
  )
  stats$mean = stats$elapsed_sec[, "mean"]
  stats$sd = stats$elapsed_sec[, "sd"]
  stats$n = stats$elapsed_sec[, "n"]
  stats$elapsed_sec = NULL
  if ("inner_repeats" %in% names(df)) {
    inner = aggregate(
      inner_repeats ~ name + label + version,
      data = df,
      FUN = function(x) unique(x)[1]
    )
    stats = merge(stats, inner, by = c("name", "label", "version"), all = TRUE)
  }
  stats
}

prev_commit = get_arg("--commit", PREV_COMMIT)
output_path = get_arg("--output")
keep_temp = has_flag("--keep-temp")

this_script = script_path()
if (is.null(this_script)) {
  stop("Unable to resolve script path; run this file with Rscript.")
}
pkg_dir = normalizePath(file.path(dirname(this_script), ".."))
repo_root = normalizePath(file.path(pkg_dir, ".."))
if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  stop("Package directory not found at: ", pkg_dir)
}

git_bin = require_tool("git")
ant_bin = require_tool("ant")
r_bin = file.path(R.home("bin"), "R")
if (!file.exists(r_bin)) {
  r_bin = require_tool("R")
}
rscript_bin = file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript_bin)) {
  rscript_bin = require_tool("Rscript")
}

old_wd = getwd()
setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

status = system2(git_bin, c("cat-file", "-e", paste0(prev_commit, "^{commit}")))
if (status != 0) {
  stop("Commit not found in local git history: ", prev_commit)
}

temp_root = tempfile("ged_compare_")
dir.create(temp_root)
prev_archive = file.path(temp_root, "prev.tar")
prev_root = file.path(temp_root, "prev_repo")
dir.create(prev_root)

cleanup_paths = c(temp_root)
if (!keep_temp) {
  on.exit(unlink(cleanup_paths, recursive = TRUE, force = TRUE), add = TRUE)
}

status = system2(git_bin, c("archive", "--format=tar", prev_commit), stdout = prev_archive)
if (status != 0) {
  stop("git archive failed for commit: ", prev_commit)
}
utils::untar(prev_archive, exdir = prev_root)

prev_pkg_dir = file.path(prev_root, "GreedyExperimentalDesign")
if (!file.exists(file.path(prev_pkg_dir, "DESCRIPTION"))) {
  stop("Previous package directory not found at: ", prev_pkg_dir)
}

lib_current = file.path(temp_root, "lib_current")
lib_prev = file.path(temp_root, "lib_prev")

cat("Rebuilding Java jar for current version...\n")
run_ant_dist(repo_root, ant_bin)
cat("Rebuilding Java jar for previous version...\n")
run_ant_dist(prev_root, ant_bin)

cat("Installing current package...\n")
install_pkg(pkg_dir, lib_current, r_bin)
cat("Installing previous package...\n")
install_pkg(prev_pkg_dir, lib_prev, r_bin)

bench_script = file.path(pkg_dir, "benchmarks", "run_version_benchmarks.R")
if (!file.exists(bench_script)) {
  stop("Benchmark runner not found at: ", bench_script)
}

current_csv = file.path(temp_root, "current_results.csv")
prev_csv = file.path(temp_root, "prev_results.csv")

cat("Running benchmarks for current version...\n")
run_bench(rscript_bin, bench_script, lib_current, current_csv, "current")
cat("Running benchmarks for previous version...\n")
run_bench(rscript_bin, bench_script, lib_prev, prev_csv, "previous")

current = read.csv(current_csv, stringsAsFactors = FALSE)
previous = read.csv(prev_csv, stringsAsFactors = FALSE)

current_stats = summarize_runs(current)
previous_stats = summarize_runs(previous)

comparison = merge(
  previous_stats,
  current_stats,
  by = "name",
  suffixes = c("_prev", "_current"),
  all = TRUE
)
comparison$ratio = comparison$mean_current / comparison$mean_prev
comparison$delta_sec = comparison$mean_current - comparison$mean_prev
comparison$se = sqrt((comparison$sd_prev^2 / comparison$n_prev) +
  (comparison$sd_current^2 / comparison$n_current))
comparison$z = ifelse(comparison$se > 0, comparison$delta_sec / comparison$se, NA_real_)
comparison$p_value = ifelse(comparison$se > 0, 2 * stats::pnorm(-abs(comparison$z)), NA_real_)
comparison$faster = comparison$ratio < 1

comparison = comparison[order(comparison$ratio), ]
print(comparison, row.names = FALSE)

if (!is.null(output_path) && output_path != "") {
  write.csv(comparison, output_path, row.names = FALSE)
  cat("Wrote comparison results to:", output_path, "\n")
}
