options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, doParallel, tidyverse, magrittr, data.table, viridis, RColorBrewer, 
               ggsci, lmtest, sandwich, ggpubr)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign



n = 100
Nw = 5000
all_w = complete_randomization_with_forced_balanced(n, Nw, form = "pos_one_min_one")

ps = c(1, 5, 10)
covariate_distributions = c("uniform", "normal", "exponential")
kernels = c("mahalanobis", "quadratic", "exponential", "gaussian")
Nsim = 25

all_res = data.table()

for (i_p in 1 : length(ps)){
  p = ps[i_p]
  for (i_c in 1 : length(covariate_distributions)){
    covariate_distribution = covariate_distributions[i_c]
    for (nsim in 1 : Nsim){
      cat("p", p, "covariate_distribution", covariate_distribution, "nsim", nsim, "\n")
      if (covariate_distribution == "uniform"){
        X = runif(n * p)
      } else if (covariate_distribution == "normal"){
        X = rnorm(n * p)
      } else if (covariate_distribution == "exponential"){
        X = rexp(n * p)
      }
      X = matrix(X, nrow = n, ncol = p)
      X = scale(X)
      
      for (i_k in 1 : length(kernels)){
        kernel = kernels[i_k]
        
        if (kernel == "mahalanobis"){
          K = X %*% solve(var(X)) %*% t(X)
        } else {
          K = matrix(NA, n, n)
          for (i in 1 : n){
            for (j in 1 : n){
              xi = X[i, , drop = FALSE]
              xj = X[j, , drop = FALSE]
              if (kernel == "quadratic"){
                K[i, j] = (1 + xi %*% t(xj) / 2)^2
              } else if (kernel == "exponential"){
                K[i, j] = exp(xi %*% t(xj))
              } else if (kernel == "gaussian"){
                K[i, j] = exp(-sum((xi - xj)^2))
              }  
            }
          }       
        }
        
        ####now run the simulation!
        obj_vals = array(NA, nrow(all_w))
        for (i_w in 1 : nrow(all_w)){
          w = all_w[i_w, , drop = FALSE]
          obj_vals[i_w] = w %*% K %*% t(w)
        }
        
        
        all_res = rbind(all_res, data.table(obj_vals = obj_vals, nsim = nsim, kernel = kernel, covariate_distribution = covariate_distribution, p = p))
      }
    }
  }
}
save(all_res, file = "kernel_large_n_bcrd_sims.RData")
all_res[, seq_var := 1 : .N, by = kernel]
all_res_ker = dcast(all_res, seq_var + nsim + p + covariate_distribution ~ kernel, value.var = c("obj_vals"))

all_corr_matrices = list()
for (p0 in ps){
  all_corr_matrices[[as.character(p0)]] = list()
  for (covariate_distribution0 in covariate_distributions){
    corr_matrix = matrix(0, nrow = length(kernels), ncol = length(kernels))
    for (nsim0 in 1 : Nsim){
      corr_matrix = corr_matrix +
        cor(as.matrix(all_res_ker[nsim == nsim0 & covariate_distribution == covariate_distribution0 & p == p0, 
                           .(exponential,  gaussian, mahalanobis, quadratic)]))      
    }
    all_corr_matrices[[as.character(p0)]][[covariate_distribution0]] = corr_matrix / Nsim
  }
}
all_corr_matrices

###viz
ggplot(all_res[nsim == 1 & p == 1 & covariate_distribution == "normal"]) + 
  geom_histogram(aes(x = log10(obj_vals), fill = kernel), alpha = 0.5, bins = 1000)
ggplot(all_res[nsim <= 25]) + 
  geom_histogram(aes(x = log10(obj_vals), 
    fill = covariate_distribution), 
    alpha = 0.5, 
    bins = 100) + 
  facet_wrap(p ~ kernel, scales = "free",
    labeller = label_wrap_gen(multi_line=FALSE)) + #labeller(p = label_both, model = label_both)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


##
p0 = 10
nsim0 = 1
covariate_distribution0 = "normal"
fig1 = ggplot(all_res_ker[nsim == nsim0 & p == p0 & covariate_distribution == covariate_distribution0]) +
  geom_point(aes(x = exponential, y = gaussian)) +
  scale_x_log10() + scale_y_log10()

fig2 = ggplot(all_res_ker[nsim == nsim0 & p == p0 & covariate_distribution == covariate_distribution0]) +
  geom_point(aes(x = mahalanobis, y = quadratic)) +
  scale_x_log10() + scale_y_log10()

fig3 = ggplot(all_res_ker[nsim == nsim0 & p == p0 & covariate_distribution == covariate_distribution0]) +
  geom_point(aes(x = mahalanobis, y = gaussian)) +
  scale_x_log10() + scale_y_log10()

fig4 = ggplot(all_res_ker[nsim == nsim0 & p == p0 & covariate_distribution == covariate_distribution0]) +
  geom_point(aes(x = mahalanobis, y = exponential)) +
  scale_x_log10() + scale_y_log10()

fig1234 = ggarrange(fig1, fig2, fig3, fig4, nrow = 2, ncol = 2)
annotate_figure(fig1234, top = paste0("nsim = ", nsim0, ", p = ", p0, ", covariate_distribution = ", covariate_distribution0))

# gd = initGreedyExperimentalDesignObject(X, max_designs = w_max, Kgram = K, objective = "kernel", diagnostics = TRUE)
# gd_res = resultsGreedySearch(gd, max_vectors = w_max)
# gd_res_ws = gd_res$ending_indicTs
# gd_res_ws[1, , drop = FALSE] %*% K %*% t(gd_res_ws[1, , drop = FALSE])
# gd_res_ws[2, , drop = FALSE] %*% K %*% t(gd_res_ws[2, , drop = FALSE])
#   
# 
# gd_res$obj_vals