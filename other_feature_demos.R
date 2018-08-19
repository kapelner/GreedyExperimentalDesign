

####Gurobi testing stuff
options(java.parameters = "-Xmx4000m")
NUM_CORES = 3
library(GreedyExperimentalDesign)
.jaddClassPath("/gurobi801/win64/lib/gurobi.jar")


#sample a covariate matrix
n = 30
p = 1
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
r = 100

#numerical optimization
gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(X, time_limit_min = 3, 
                                                                num_cores = NUM_CORES)
indicTs = resultsGurobiNumericalOptimizeExperimentalDesign(gnoed)$indicTs
indicTs
compute_objective_val(X, indicTs[,1], objective = "mahal_dist")

random_indices = sample(1 : n)
X_randomized = X[random_indices, , drop = FALSE]
gnoed = initGurobiNumericalOptimizationExperimentalDesignObject(X_randomized, time_limit_min = 0.5, num_cores = NUM_CORES)
indicT = resultsGurobiNumericalOptimizeExperimentalDesign(gnoed)$indicTs[,1]
compute_objective_val(X_randomized, indicT, objective = "mahal_dist")
compute_objective_val(X, indicT[order(random_indices)], objective = "mahal_dist")


#Gurobi: eigenvalue estimation of var-cov matrix
indicTs = gurobi_multiple_designs(X, r, time_limit_min = 0.5, num_cores = NUM_CORES)
indicTs = 2 * (indicTs - 0.5) #make it -1, 1
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]

#greedy pair switching: eigenvalue estimation of var-cov matrix
options(java.parameters = "-Xmx4000m")
NUM_CORES = 3
library(GreedyExperimentalDesign)

rd = initGreedyExperimentalDesignObject(X, r, wait = TRUE, objective = "mahal_dist", num_cores = 4)
res = resultsGreedySearch(rd, max_vectors = r)
indicTs = t(res$ending_indicTs)
for (i in 1 : nrow(indicTs)){
  indicTs[i, ] = 2 * (indicTs[i, ] - 0.5)
}
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]
# res$obj_vals


#CR: eigenvalue estimation of var-cov matrix
indicTs = complete_randomization_with_forced_balanced(n, r)
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]

indicTs = complete_randomization(n, r)
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]





#### Kernel stuff
options(java.parameters = "-Xmx4000m")
NUM_CORES = 3
library(GreedyExperimentalDesign)

n = 20
p = 3
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
Kgram = compute_gram_matrix(X, "rbf", 1)

#check it works for greedy
rd = initGreedyExperimentalDesignObject(wait = TRUE, objective = "kernel", num_cores = 1, Kgram = Kgram)
res = resultsGreedySearch(rd, max_vectors = r)
indicTs = unique(t(res$ending_indicTs))
for (i in 1 : nrow(indicTs)){
  indicTs[i, ] = 2 * (indicTs[i, ] - 0.5)
}
var_cov_matrix_est = var(indicTs)
abs(var_cov_matrix_est[1:10, 1:10])
eigen(var_cov_matrix_est)$values[1]

#check it works for optimal
id = initOptimalExperimentalDesignObject(Kgram = Kgram, objective = "kernel", wait = TRUE)
indicT_kernel = resultsOptimalSearch(id)$indicT
id = initOptimalExperimentalDesignObject(X = X, objective = "mahal_dist", wait = TRUE)
indicT_mahal = resultsOptimalSearch(id)$indicT
compute_objective_val(X, indicT_kernel, objective = "mahal_dist")
compute_objective_val(X, indicT_mahal, objective = "mahal_dist")

#check it works for rerandomization
rd = initRerandomizationExperimentalDesignObject(Kgram = Kgram, objective = "kernel", wait = TRUE)
ending_indicTs = resultsRerandomizationSearch(rd, include_assignments = TRUE)$ending_indicTs
dim(ending_indicTs)
rd = initRerandomizationExperimentalDesignObject(Kgram = Kgram, objective = "kernel", wait = TRUE, obj_val_cutoff_to_include = 18)
head(sort(resultsRerandomizationSearch(rd)$obj_vals), 200)
ending_indicTs = resultsRerandomizationSearch(rd, include_assignments = TRUE)$ending_indicTs
dim(ending_indicTs)


#check it works for Gurobi
options(java.parameters = "-Xmx4000m")
NUM_CORES = 3
library(GreedyExperimentalDesign)
.jaddClassPath("/gurobi801/win64/lib/gurobi.jar")

n = 40
p = 3
X = generate_stdzied_design_matrix(n = n, p = p, covariate_gen = rnorm)
Kgram = compute_gram_matrix(X, "rbf", 1)


rd = initGurobiNumericalOptimizationExperimentalDesignObject(Kgram = Kgram, objective = "kernel", 
                                                             time_limit_min = 3, 
                                                             num_cores = NUM_CORES)
indicT_kernel_rbf = resultsGurobiNumericalOptimizeExperimentalDesign(rd)$indicTs[, 1]

Kgram = compute_gram_matrix(X, "poly", c(1, 1, 1)) #linear

rd = initGurobiNumericalOptimizationExperimentalDesignObject(Kgram = Kgram, objective = "kernel", 
                                                             time_limit_min = 3, 
                                                             num_cores = NUM_CORES)
indicT_kernel_linear = resultsGurobiNumericalOptimizeExperimentalDesign(rd)$indicTs[, 1]


rd = initGurobiNumericalOptimizationExperimentalDesignObject(X, objective = "mahal_dist", 
                                                             time_limit_min = 3, 
                                                             num_cores = NUM_CORES)
indicT_mahal = resultsGurobiNumericalOptimizeExperimentalDesign(rd)$indicTs[, 1]

indicT_kernel_rbf
indicT_kernel_linear
indicT_mahal
compute_objective_val(X, indicT_kernel_rbf, objective = "mahal_dist")
compute_objective_val(X, indicT_kernel_linear, objective = "mahal_dist")
compute_objective_val(X, indicT_mahal, objective = "mahal_dist")








