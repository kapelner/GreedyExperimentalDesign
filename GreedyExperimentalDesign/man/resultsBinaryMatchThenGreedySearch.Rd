% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_match_followed_by_greedy_search.R
\name{resultsBinaryMatchThenGreedySearch}
\alias{resultsBinaryMatchThenGreedySearch}
\title{Returns unique allocation vectors that are binary matched}
\usage{
resultsBinaryMatchThenGreedySearch(
  obj,
  num_vectors = NULL,
  compute_obj_vals = FALSE,
  form = "zero_one"
)
}
\arguments{
\item{obj}{The \code{binary_then_greedy_experimental_design} object where the pairs are computed.}

\item{num_vectors}{How many random allocation vectors you wish to return. The default is \code{NULL} indicating you want all of them.}

\item{compute_obj_vals}{Should we compute all the objective values for each allocation? Default is \code{FALSE}.}

\item{form}{Which form should it be in? The default is \code{one_zero} for 1/0's or \code{pos_one_min_one} for +1/-1's.}
}
\description{
Returns unique allocation vectors that are binary matched
}
\author{
Adam Kapelner
}
