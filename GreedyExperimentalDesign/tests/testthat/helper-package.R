options(java.parameters = "-Xmx10g")

library(testthat)

MAX_FAILS = 3
testthat::set_max_fails(MAX_FAILS)
options(testthat.progress.max_fails = MAX_FAILS)
Sys.setenv(TESTTHAT_MAX_FAILS = as.character(MAX_FAILS))
Sys.setenv(TESTTHAT_PARALLEL = "false")

ged_failure_state = new.env(parent = emptyenv())
ged_failure_state$count = 0L
ged_failure_state$max_fail = MAX_FAILS
ged_failure_state$aborting = FALSE

resolve_testthat_dir = function() {
  wd = getwd()
  if (dir.exists(file.path(wd, "tests", "testthat"))) {
    return(file.path(wd, "tests", "testthat"))
  }
  if (dir.exists(file.path(wd, "GreedyExperimentalDesign", "tests", "testthat"))) {
    return(file.path(wd, "GreedyExperimentalDesign", "tests", "testthat"))
  }
  if (basename(wd) == "testthat" && basename(dirname(wd)) == "tests") {
    return(wd)
  }
  wd
}

resolve_package_root = function() {
  testthat_dir = resolve_testthat_dir()
  candidate = normalizePath(file.path(testthat_dir, "..", ".."), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(candidate, "DESCRIPTION"))) {
    return(candidate)
  }
  candidate
}

package_root = resolve_package_root()
local_lib = file.path(package_root, ".Rlib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

add_java_classpath = function(root_dir) {
  jar_path = file.path(root_dir, "inst", "java", "GreedyExperimentalDesign.jar")
  if (file.exists(jar_path)) {
    rJava::.jaddClassPath(jar_path)
  }
  invisible(NULL)
}

use_pkgload = requireNamespace("pkgload", quietly = TRUE) &&
  file.exists(file.path(package_root, "DESCRIPTION")) &&
  dir.exists(file.path(package_root, "R"))
if (use_pkgload) {
  pkgload::load_all(package_root, export_all = FALSE, quiet = TRUE)
  add_java_classpath(package_root)
} else {
  library(GreedyExperimentalDesign)
}

ged_failure_state$log_path = file.path(resolve_testthat_dir(), "ged-failures.log")
dir.create(dirname(ged_failure_state$log_path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(ged_failure_state$log_path)) {
  file.remove(ged_failure_state$log_path)
}
cat(
  paste0("---- ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ----\n"),
  "Test run started.\n\n",
  file = ged_failure_state$log_path,
  append = TRUE
)

append_failure_log = function(exp) {
  if (is.null(ged_failure_state$log_path) || !nzchar(ged_failure_state$log_path)) {
    ged_failure_state$log_path = file.path(tempdir(), "ged-failures.log")
  }
  header = paste0(
    "---- ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " ----\n"
  )
  body = testthat:::issue_summary(exp)
  cat(
    header,
    body,
    "\n\n",
    file = ged_failure_state$log_path,
    append = TRUE
  )
}

emit_failure_console = function(exp) {
  header = paste0(
    "---- ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " ----\n"
  )
  body = testthat:::issue_summary(exp)
  cat(header, body, "\n\n", file = stderr())
  flush.console()
}

record_immediate_failure = function(exp) {
  ged_failure_state$count = ged_failure_state$count + 1L
  append_failure_log(exp)
  emit_failure_console(exp)
  if (ged_failure_state$count >= ged_failure_state$max_fail) {
    ged_failure_state$aborting = TRUE
    stop("Too many failures; aborting remaining tests.")
  }
  invisible(NULL)
}

LoggingReporter = R6::R6Class(
  "LoggingReporter",
  inherit = testthat::Reporter,
  public = list(
    add_result = function(context, test, result) {
      if (!testthat:::expectation_success(result)) {
        record_immediate_failure(result)
      }
    }
  )
)

install_reporter_hook = function() {
  if (!exists("ged_orig_test_files_reporter", envir = .GlobalEnv, inherits = FALSE)) {
    assign("ged_orig_test_files_reporter", testthat:::test_files_reporter, envir = .GlobalEnv)
    utils::assignInNamespace(
      "test_files_reporter",
      function(reporter, mode = c("serial", "parallel"), desc = NULL, frame = rlang::caller_env()) {
        res = get("ged_orig_test_files_reporter", envir = .GlobalEnv)(
          reporter,
          mode = mode,
          desc = desc,
          frame = frame
        )
        res$multi = testthat::MultiReporter$new(reporters = c(list(LoggingReporter$new()), res$multi$reporters))
        res
      },
      ns = "testthat"
    )
  }
  invisible(TRUE)
}

patch_recover2 = function() {
  if (!exists("ged_orig_recover2", envir = .GlobalEnv, inherits = FALSE)) {
    assign("ged_orig_recover2", testthat:::recover2, envir = .GlobalEnv)
    utils::assignInNamespace(
      "recover2",
      function(...) {
        if (!interactive()) {
          return(invisible(FALSE))
        }
        get("ged_orig_recover2", envir = .GlobalEnv)(...)
      },
      ns = "testthat"
    )
  }
  invisible(TRUE)
}

abort_if_too_many_failures = function() {
  if (isTRUE(ged_failure_state$aborting)) {
    stop("Too many failures; aborting remaining tests.")
  }
  invisible(NULL)
}

with_immediate_failures = function(expr) {
  abort_if_too_many_failures()
  force(expr)
}

skip_on_cmd_check = function() {
  is_check = nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) ||
    nzchar(Sys.getenv("R_CMD_CHECK")) ||
    nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_")) ||
    nzchar(Sys.getenv("_R_CHECK_CRAN_INCOMING_REMOTE_"))
  if (is_check) {
    testthat::skip("Skipped during R CMD check.")
  }
  invisible(NULL)
}

install_reporter_hook()
patch_recover2()
