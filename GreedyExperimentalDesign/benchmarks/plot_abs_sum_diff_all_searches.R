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
if (is.null(this_script)) {
  stop("Unable to resolve script path; run this file with Rscript.")
}

script_dir = dirname(this_script)
target = file.path(script_dir, "plot_obj_value_comparison_all_searches.R")
if (!file.exists(target)) {
  stop("Expected script not found: ", target)
}

source(target)
