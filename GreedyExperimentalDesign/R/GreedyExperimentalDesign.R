#' A tool to find many types of a priori experimental designs
#'
#' @name 		GreedyExperimentalDesign
#' @title 		Greedy Experimental Design Search
#' @author 		Adam Kapelner \email{kapelner@@qc.cuny.edu}
#' @references 	Kapelner, A
#' @keywords 	design optimize
#' @import      rJava stats graphics grDevices checkmate
#' @importFrom 	nbpMatching distancematrix nonbimatch
#' @importFrom  Rcpp sourceCpp evalCpp
#' @importFrom  stringr str_detect
#' @importFrom  stringi stri_replace_first_regex
#' @importFrom  rlist list.cbind
#' @importFrom  checkmate assertClass assertCount assertNumeric assertChoice assertLogical assertTRUE assertMatrix vname
#' @useDynLib 	GreedyExperimentalDesign
##### Run "library(roxygen2); roxygenise("GreedyExperimentalDesign", clean = TRUE)" to regenerate all Rd files and NAMESPACE and DESCRIPTION file
##### but make sure you are in the root directory of the project. Make sure to add these two to namespace afterwards:
NULL