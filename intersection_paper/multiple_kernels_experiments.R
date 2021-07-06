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
p = 5
Nw = 2000

covariate_distribution = "uniform"

if (covariate_distribution == "uniform"){
  X = runif(n * p)
} else if (covariate_distribution == "normal"){
  X = rnorm(n * p)
} else if (covariate_distribution == "exponential"){
  X = rexp(n * p)
}
X = matrix(X, nrow = n, ncol = p)
X = scale(X)



gmked = initGreedyMultipleKernelExperimentalDesignObject(X, 
   max_designs = Nw, kernel_pre_num_designs = Nw, seed = 1, maximum_gain_scaling = 1,
   kernel_names = kernel_names, diagnostics = TRUE, wait = TRUE, num_cores = nC)
# startSearch(gmked)
gmked_res = resultsMultipleKernelGreedySearch(gmked, max_vectors = Nw, form = "pos_one_min_one")
# gmked_res$obj_vals %>% mean
# gmked_res$num_iters
# gmked_res$kernel_obj_vals_by_iter
# gmked_res$obj_vals
# lapply(gmked_res$kernel_obj_vals_by_iter, rowMeans)


ggplot(data.frame(obj_vals = gmked_res$obj_vals)) + 
  geom_histogram(aes(x = obj_vals), bins = 100)

ggplot_obj = ggplot(do.call("rbind", gmked$all_univariate_kernel_data))
ggplot_obj +
  geom_histogram(aes(x = log10objvalsi), fill = "red", bins = 50) +
  facet_wrap(kernel ~ ., scales = "free")

ggplot_obj +
  geom_histogram(aes(x = log10objvalsf), fill = "red", bins = 50) +
  facet_wrap(kernel ~ ., scales = "free")

ggplot_obj +
  geom_histogram(aes(x = log10_i_over_f), fill = "red", bins = 50) +
  facet_wrap(kernel ~ ., scales = "free")

ggplot_obj +
  geom_histogram(aes(x = pct_red_max), fill = "red", bins = 50) +
  facet_wrap(kernel ~ ., scales = "free")




Kgrams = list()
all_res = data.table()
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
        #http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#inverse_multiquadric
        if (kernel == "poly_2"){
          Kgrams[[i_k]][i, j] = (1 + xi %*% t(xj) / 2)^2
        } else if (kernel == "exponential"){
          Kgrams[[i_k]][i, j] = exp(xi %*% t(xj))
        } else if (kernel == "gaussian"){
          Kgrams[[i_k]][i, j] = exp(-sum((xi - xj)^2))
        } else if (kernel == "laplacian"){
          Kgrams[[i_k]][i, j] = exp(-sqrt(sum((xi - xj)^2)))
        } else if (kernel == "inv_mult_quad"){
          Kgrams[[i_k]][i, j] = 1 / sqrt(sum((xi - xj)^2) + 1)
        }
      }
    }       
  }
  
  cat("kernel:", kernel, "\n")
  gd = initGreedyExperimentalDesignObject(X, 
          max_designs = Nw, Kgram = Kgrams[[i_k]], objective = "kernel",
          diagnostics = TRUE, wait = TRUE, seed = 1, num_cores = nC)
  gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
  starting_obj_vals = array(NA, Nw)
  for (i in 1 : Nw){
    starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
  }
  all_res = rbind(all_res, data.table(strategy = "solo", kernel = kernel, s = starting_obj_vals, e = gd_res$obj_vals_unordered))
  
  ending_obj_vals = array(NA, Nw)
  for (i in 1 : Nw){
    starting_obj_vals[i] = gmked_res$starting_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gmked_res$starting_indicTs[i, , drop = FALSE])
    ending_obj_vals[i] = gmked_res$ending_indicTs[i, , drop = FALSE] %*% Kgrams[[i_k]] %*% t(gmked_res$ending_indicTs[i, , drop = FALSE])
  }
  gd_res$obj_vals_unordered
  all_res = rbind(all_res, data.table(strategy = "multi_pct", kernel = kernel, s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_sum", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_sum_norm", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_prod_norm", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_sum_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
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
  all_res = rbind(all_res, data.table(strategy = "multi_sum_plus_prod", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
}


################# INVERSE VARIANCE WEIGHTING!!!
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
  all_res = rbind(all_res, data.table(strategy = "multi_weighted_sum", kernel = kernel_names[i_k], s = starting_obj_vals, e = ending_obj_vals))
}

all_res %<>%
  mutate(log_s = log10(s)) %>%
  mutate(log_e = log10(e)) %>%
  mutate(log_s_e = log10(s / e)) %>%
  group_by(kernel) %>%
  mutate(pct_reduction_lacking = 1 - log_s_e / (max(log_s_e) * 1.1)) %>%
  arrange(desc(log_s_e)) %>%
  ungroup
all_res

all_res_wide = all_res %>% pivot_longer(cols = c(s, e, log_s, log_e, log_s_e, pct_reduction_lacking))

ggplot(all_res_wide %>% filter(name %in% c("log_s", "log_e"))) +
  geom_histogram(aes(x = value, fill = strategy), alpha = 0.5, bins = 100) + 
  facet_wrap(kernel ~ ., scales = "free_y") + xlim(-13, 5)

ggplot(all_res) +
  geom_histogram(aes(x = log_s_e, fill = strategy), alpha = 0.5, bins = 100) +
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

ggplot(all_res) +
  geom_density(aes(x = pct_reduction_lacking, fill = strategy), alpha = 0.5) +
  facet_grid(kernel ~ ., scales = "free_y")

all_res_summary = all_res %>%
  group_by(strategy, kernel) %>%
  summarize(avg_log_s_e = mean(log_s_e), upper_10 = quantile(log_s_e, .9)) %>% 
  arrange(kernel, desc(avg_log_s_e)) %>%
  ungroup %>%
  group_by(kernel) %>%
  mutate(rank = 1 : n()) %>%
  mutate(pct_worse_than_best = (first(avg_log_s_e) - avg_log_s_e) / first(avg_log_s_e) * 100)

all_res_summary %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  as.data.frame
with(all_res_summary, table(strategy, rank))

all_res_summary %>%
  group_by(strategy) %>%
  summarize(avg_pct_worse_than_best = mean(pct_worse_than_best)) %>%
  arrange(avg_pct_worse_than_best) %>%
  mutate_if(is.numeric, round, digits=3) %>%
  as.data.frame


#




























ggplot(data.frame(s = starting_obj_vals, e = gd_res$obj_vals_unordered)) + 
  geom_point(aes(x = s, y = log10(e)))

gd = initGreedyExperimentalDesignObject(X, 
      max_designs = Nw, Kgram = K, objective = "kernel",
      diagnostics = TRUE, wait = TRUE, num_cores = nC)
gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
gd_res_ws = gd_res$ending_indicTs
starting_obj_vals = array(NA, Nw)
ending_obj_vals = array(NA, Nw)
for (i in 1 : Nw){
  starting_obj_vals[i] = gd_res$starting_indicTs[i, , drop = FALSE] %*% K %*% t(gd_res$starting_indicTs[i, , drop = FALSE])
  ending_obj_vals[i] = gd_res_ws[i, , drop = FALSE] %*% K %*% t(gd_res_ws[i, , drop = FALSE])
}
obj_vals_0 = starting_obj_vals[order(ending_obj_vals)]
obj_vals_f = as.numeric(ending_obj_vals[order(ending_obj_vals)])
gd_res$obj_vals
all.equal(obj_vals_f, gd_res$obj_vals)

hist(obj_vals_0, breaks = 50)
hist(obj_vals_f, breaks = 50)
hist(log10(obj_vals_0 / obj_vals_f), breaks = 50)

# gd_res$obj_val_by_iters
# gd_res$num_iters

# obj_val_by_iters = list()
# for (i in 1 : Nw){
#   switches = gd_res$switches[[i]] + 1
#   obj_val_by_iters[[i]] = array(NA, nrow(switches))
#   w = gd_res$starting_indicTs[i, , drop = FALSE]
#   obj_val_by_iters[[i]][1] = w %*% K %*% t(w)
#   for (s in 1 : nrow(switches)){
#     temp_wl = w[switches[s, 1]]
#     w[switches[s, 1]] = w[switches[s, 2]]
#     w[switches[s, 2]] = -w[switches[s, 1]]
#     assert(sum(w) == 0)
#     obj_val_by_iters[[i]][s + 1] = w %*% K %*% t(w) #eval_quad_form(K, w)
#   }
# }
# obj_val_by_iters
# gd_res$obj_val_by_iters


# eval_quad_form = function(A, x){
#   tot = 0
#   for (i in 1 : nrow(A)){
#     for (j in 1 : ncol(A)){
#       tot = tot + A[i, j] * x[i] * x[j]
#     }
#   }
#   tot
# }










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
              #https://people.eecs.berkeley.edu/~jordan/kernels/
              #http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#inverse_multiquadric
              if (kernel == "quadratic"){
                K[i, j] = (1 + xi %*% t(xj) / 2)^2
              } else if (kernel == "exponential"){
                K[i, j] = exp(xi %*% t(xj))
              } else if (kernel == "gaussian"){
                K[i, j] = exp(-sum((xi - xj)^2))
              } else if (kernel == "laplacian"){
                K[i, j] = exp(-sqrt(sum((xi - xj)^2)))
              } else if (kernel == "inv_mult_quad"){
                K[i, j] = 1 / sqrt(sum((xi - xj)^2) + 1)
              }
            }
          }       
        }
        
        # Nw = 1000
        # microbenchmark::microbenchmark(
        #   init = {initGreedyExperimentalDesignObject(X, max_designs = Nw, Kgram = K, objective = "kernel", diagnostics = FALSE, wait = TRUE, num_cores = nC)},
        #   times = 2
        # )
        # for (i in 1 : 100){
        #   gd = initGreedyExperimentalDesignObject(X, max_designs = Nw, Kgram = K, objective = "kernel", diagnostics = FALSE, wait = TRUE, num_cores = nC)
        # }
        
        
        gd = initGreedyExperimentalDesignObject(X, max_designs = Nw, Kgram = K, objective = "kernel", diagnostics = TRUE, wait = TRUE, num_cores = nC)
        gd_res = resultsGreedySearch(gd, max_vectors = Nw, form = "pos_one_min_one")
        gd_res_ws = gd_res$ending_indicTs
        gd_res_ws[1, , drop = FALSE] %*% K %*% t(gd_res_ws[1, , drop = FALSE])
        gd_res_ws[2, , drop = FALSE] %*% K %*% t(gd_res_ws[2, , drop = FALSE])


        gd_res$obj_vals
        
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

