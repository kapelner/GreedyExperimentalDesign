options(java.parameters = c("-Xmx20000m"))
pacman::p_load(GreedyExperimentalDesign, tidyverse, magrittr, data.table, ggpubr, microbenchmark, checkmate)
#R CMD INSTALL -l ~/Documents/R/win-library/3.6/ GreedyExperimentalDesign

nC = 5
covariate_distributions = c("uniform", "normal", "exponential")
kernel_names = c("mahalanobis", "poly_2", "exponential", "gaussian", "laplacian", "inv_mult_quad")
m = length(kernel_names)


ns = c(32, 100, 132, 300)
ps = c(1, 2, 5, 10)
n = 100
Nw = 2000

all_res = data.table()
for (covariate_distribution in covariate_distributions){
  for (p in ps){
    cat("covariate_distribution:", covariate_distribution, " p:", p, "\n")
    
    if (covariate_distribution == "uniform"){
      X = runif(n * p)
    } else if (covariate_distribution == "normal"){
      X = rnorm(n * p)
    } else if (covariate_distribution == "exponential"){
      X = rexp(n * p)
    }
    X = matrix(X, nrow = n, ncol = p)
    X = scale(X)
    
    # gmked = initGreedyMultipleKernelExperimentalDesignObject(X, 
    #    max_designs = Nw, kernel_pre_num_designs = Nw, seed = 1, maximum_gain_scaling = 1,
    #    kernel_names = kernel_names, diagnostics = TRUE, wait = TRUE, num_cores = nC)
    # # startSearch(gmked)
    # gmked_res = resultsMultipleKernelGreedySearch(gmked, max_vectors = Nw, form = "pos_one_min_one")
    # # gmked_res$obj_vals %>% mean
    # # gmked_res$num_iters
    # # gmked_res$kernel_obj_vals_by_iter
    # # gmked_res$obj_vals
    # # lapply(gmked_res$kernel_obj_vals_by_iter, rowMeans)
    # 
    # 
    # ggplot(data.frame(obj_vals = gmked_res$obj_vals)) + 
    #   geom_histogram(aes(x = obj_vals), bins = 100)
    # 
    # ggplot_obj = ggplot(do.call("rbind", gmked$all_univariate_kernel_data))
    # ggplot_obj +
    #   geom_histogram(aes(x = log10objvalsi), fill = "red", bins = 50) +
    #   facet_wrap(kernel ~ ., scales = "free")
    # 
    # ggplot_obj +
    #   geom_histogram(aes(x = log10objvalsf), fill = "red", bins = 50) +
    #   facet_wrap(kernel ~ ., scales = "free")
    # 
    # ggplot_obj +
    #   geom_histogram(aes(x = log10_i_over_f), fill = "red", bins = 50) +
    #   facet_wrap(kernel ~ ., scales = "free")
    # 
    # ggplot_obj +
    #   geom_histogram(aes(x = pct_red_max), fill = "red", bins = 50) +
    #   facet_wrap(kernel ~ ., scales = "free")
    
    Kgrams = list()
    for (i_k in 1 : m){
      kernel = kernel_names[i_k]
      
      if (kernel == "mahalanobis"){
        Kgrams[[i_k]] = X %*% solve(var(X)) %*% t(X)
      } else {
        Kgrams[[i_k]] = matrix(NA, n, n)
        for (i in 1 : n){
          for (j in 1 : n){
            xi = X[i, , drop = FALSE]
            xj = X[j, , drop = FALSE]
            #https://people.eecs.berkeley.edu/~jordan/kernels/0521813972c01_p03-24.pdf
            if (kernel == "poly_2"){
              Kgrams[[i_k]][i, j] = (1 + xi %*% t(xj) / 2)^2
            } else if (kernel == "exponential"){
              Kgrams[[i_k]][i, j] = exp(xi %*% t(xj))
            } else if (kernel == "gaussian"){
              Kgrams[[i_k]][i, j] = exp(-sum((xi - xj)^2))
            } else if (kernel == "laplacian"){
              Kgrams[[i_k]][i, j] = exp(-sqrt(sum((xi - xj)^2)))
            } else if (kernel == "inv_mult_quad"){
              #http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#inverse_multiquadric
              Kgrams[[i_k]][i, j] = 1 / sqrt(sum((xi - xj)^2) + 1)
            }
          }
        }       
      }
      
      gd = initGreedyExperimentalDesignObject(X, 
              max_designs = Nw, Kgram = Kgrams[[i_k]], objective = "kernel",
              diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
      gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
      starting_obj_vals = array(NA, Nw)
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "solo", kernel = kernel, s = starting_obj_vals, e = ending_obj_vals))
      
      # for (i in 1 : Nw){
      #   starting_obj_vals[i] = gmked_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gmked_res$starting_indicTs[i, , drop = FALSE])
      #   ending_obj_vals[i] = gmked_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gmked_res$ending_indicTs[i, , drop = FALSE])
      # }
      # gd_res$obj_vals_unordered
      # all_res = rbind(all_res, data.table(strategy = "multi_pct", kernel = kernel, s = starting_obj_vals, e = ending_obj_vals))
    }
    
    #sum
    Ksum = matrix(0, n, n)
    for (i_k in 1 : m){
      Ksum = Ksum + Kgrams[[i_k]]
    }
    gd = initGreedyExperimentalDesignObject(X, 
            max_designs = Nw, Kgram = Ksum, objective = "kernel",
            diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      gd_res$obj_vals_unordered
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_sum", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    #normsum
    
    Knormsum = matrix(0, n, n)
    for (i_k in 1 : m){
      Knormsum = Knormsum + Kgrams[[i_k]] / sum(Kgrams[[i_k]])
    }
    gd = initGreedyExperimentalDesignObject(X, 
          max_designs = Nw, Kgram = Knormsum, objective = "kernel",
          diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_sum_norm", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    #prod
    
    Kprod = matrix(1, n, n)
    for (i_k in 1 : m){
      Kprod = Kprod * Kgrams[[i_k]]
    }
    gd = initGreedyExperimentalDesignObject(X, 
                                            max_designs = Nw, Kgram = Kprod, objective = "kernel",
                                            diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      gd_res$obj_vals_unordered
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    #normprod
    
    Knormprod = matrix(1, n, n)
    for (i_k in 1 : m){
      Knormprod = Knormprod * Kgrams[[i_k]] / sum(Kgrams[[i_k]])
    }
    gd = initGreedyExperimentalDesignObject(X, 
                                            max_designs = Nw, Kgram = Knormprod, objective = "kernel",
                                            diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      gd_res$obj_vals_unordered
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_prod_norm", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    
    #sum plus prod
    
    Ksumprod = matrix(0, n, n)
    for (i_k in 1 : m){
      Ksumprod = Ksumprod + Kgrams[[i_k]] * Kprod
    }
    gd = initGreedyExperimentalDesignObject(X, 
                                            max_designs = Nw, Kgram = Ksumprod, objective = "kernel",
                                            diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      gd_res$obj_vals_unordered
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_sum_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    #sum times prod
    gd = initGreedyExperimentalDesignObject(X, 
                                            max_designs = Nw, Kgram = Ksum + Kprod, objective = "kernel",
                                            diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      gd_res$obj_vals_unordered
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_sum_plus_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    
    
    ################# INVERSE SD WEIGHTING!!!
    sd_by_kernel = array(NA, m)
    for (i_k in 1 : m){
      sd_by_kernel[i_k] = all_res %>%
        filter(strategy == "solo" & kernel == kernel_names[i_k]) %>%
        mutate(log_s_e = log10(s / e)) %>%
        pull(log_s_e) %>%
        sd
    }
    sd_by_kernel = sd_by_kernel / sum(sd_by_kernel)
    
    
    Kwsum = matrix(0, n, n)
    for (i_k in 1 : m){
      Kwsum = Kwsum + Kgrams[[i_k]] * sd_by_kernel[i_k]
    }
    
    gd = initGreedyExperimentalDesignObject(X, 
          max_designs = Nw, Kgram = Kwsum, objective = "kernel",
          diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
    gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
    
    for (i_k in 1 : m){
      ending_obj_vals = array(NA, Nw)
      for (i in 1 : Nw){
        starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
        ending_obj_vals[i] = gd_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$ending_indicTs[i, , drop = FALSE])
      }
      all_res = rbind(all_res, data.table(covariate_distribution = covariate_distribution, p = p, strategy = "multi_weighted_sum", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
    }
    save(all_res, file = "mult_kern_strategies_covs_ps.RData")
  }
}

#objective functions cannot be negative so it must be numerical error
all_res %<>% filter(e > 0)

all_res %<>%
  mutate(log_s = log10(s)) %>%
  mutate(log_e = log10(e)) %>%
  mutate(log_s_e = log10(s / e)) %>%
  group_by(covariate_distribution, p) %>%
  mutate(avg_solo_performance = log_s_e * (strategy == "solo") / sum(strategy == "solo")) %>%
  group_by(covariate_distribution, p, kernel) %>%
  mutate(pct_reduction_lacking = 1 - log_s_e / avg_solo_performance) %>%
  arrange(covariate_distribution, p, kernel, desc(log_s_e)) %>%
  ungroup
all_res

all_res_wide = all_res %>% pivot_longer(cols = c(s, e, log_s, log_e, log_s_e, pct_reduction_lacking))

# ggplot(all_res_wide %>% filter(name %in% c("log_s", "log_e"))) +
#   geom_histogram(aes(x = value, fill = strategy), alpha = 0.5, bins = 100) + 
#   facet_wrap(kernel ~ p + covariate_distribution, scales = "free_y") + xlim(-13, 5)

# ggplot(all_res %>% filter(p == 1 & covariate_distribution == "normal")) +
#   geom_histogram(aes(x = log_s_e, fill = strategy), alpha = 0.5, bins = 100) +
#   facet_grid(kernel ~ ., scales = "free_y")

ggplot(all_res %>% filter(p == 1 & covariate_distribution == "uniform")) +
  geom_density(aes(x = log_s_e, fill = strategy), alpha = 0.5) +
  facet_wrap(kernel ~ ., scales = "free", ncol = 1) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

all_res %>% filter(p == 1 & covariate_distribution == "normal" & log_s_e < 0)


all_res %>% filter(covariate_distribution == "normal" & log_s_e < 0) %>%
  select(strategy, kernel, p) %>%
   table


ggplot(all_res %>% filter(p == 1 & covariate_distribution == "uniform" &
  strategy %in% c("solo", "multi_prod"))) + #, "multi_sum"
  geom_density(aes(x = log_s_e, fill = strategy), alpha = 0.5) +
  facet_wrap(kernel ~ ., scales = "free", ncol = 1) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggplot(all_res %>% filter(p == 1 & covariate_distribution == "uniform")) +
  geom_density(aes(x = log_s_e, fill = strategy), alpha = 0.5) +
  facet_grid(kernel ~ ., scales = "free_y")

ggplot(all_res) +
  geom_histogram(aes(x = pct_reduction_lacking, fill = strategy), bins = 100, alpha = 0.5) +
  facet_grid(kernel ~ ., scales = "free_y")


ggplot(all_res_wide %>% filter(name %in% c("log_s", "log_e"))) +
  geom_density(aes(x = value, fill = strategy), alpha = 0.5) + 
  facet_wrap(kernel ~ ., scales = "free_y")

ggplot(all_res %>% filter(strategy %in% c("solo", "multi_weighted_sum"))) +
  geom_density(aes(x = log_s_e, fill = strategy), alpha = 0.5) +
  facet_wrap(kernel ~ ., scales = "free")

ggplot(all_res %>% filter(p == 1 & covariate_distribution == "exponential")) +
  geom_density(aes(x = pct_reduction_lacking, fill = strategy), alpha = 0.5) +
  facet_grid(kernel ~ ., scales = "free_y")

all_res_summary = all_res %>%
  group_by(strategy, kernel, p, covariate_distribution) %>%
  summarize(avg_log_s_e = mean(log_s_e), avg_solo_log_s_e = mean(avg_solo_performance), upper_10 = quantile(log_s_e, .9)) %>% 
  arrange(kernel, desc(avg_log_s_e)) %>%
  ungroup %>%
  group_by(kernel, p, covariate_distribution) %>%
  mutate(pct_worse_than_solo = (avg_solo_log_s_e - avg_log_s_e) / avg_solo_log_s_e * 100) %>%
  mutate(rank = 1 : n())

all_res_summary %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  as.data.frame
with(all_res_summary, table(strategy, rank))


all_res_summary %>% 
  group_by(strategy) %>%
  summarize(avg_pct_worse_than_solo = mean(pct_worse_than_solo)) %>%
  arrange(avg_pct_worse_than_solo)

all_res_summary %>%
  group_by(covariate_distribution, p, strategy) %>%
  summarize(avg_pct_worse_than_solo = mean(pct_worse_than_solo)) %>%
  arrange(covariate_distribution, p, avg_pct_worse_than_solo) %>%
  mutate(rank = 1 : n()) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  as.data.frame
