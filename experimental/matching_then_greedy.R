options(java.parameters = c("-Xmx4000m")) #this means 4GB of RAM for YARF
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

X = data.matrix(MASS::Pima.tr[, 1:7])

bmeds = binaryMatchExperimentalDesignSearch(X)
bmeds$indices_pairs

bmfged = binaryMatchFollowedByGreedyExperimentalDesignSearch(X)
bmfged$binary_match_design
bmfged$greedy_design
gs_res = resultsGreedySearch(bmfged$greedy_design)
bmfged_res = resultsBinaryMatchThenGreedySearch(bmfged)
