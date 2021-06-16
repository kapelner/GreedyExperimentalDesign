
nX = 1
nEPS = 1
nR = 1e5

source("common_sim_params.R")

exp_settings = data.frame(matrix(NA, nrow = length(ps) * nX * nEPS, ncol = 3))
names(exp_settings) = c("p", "nx", "neps")
exp_settings$neps = rep(1 : nEPS, times = length(ps) * nX)
exp_settings$nx   = rep(1 : nX,   each  = nEPS)
exp_settings$p    = rep(ps,       each  = nEPS * nX)

filename = "sec_3"

source("common_sims.R")



# get ideal error which is...
Nsim_ideal = 10000000
bcrds = complete_randomization_with_forced_balanced(n, Nsim_ideal, form = "pos_one_min_one")
betaThat_ideals = array(NA, Nsim_ideal)
for (nsim_ideal in 1 : Nsim_ideal){
  w = bcrds[nsim_ideal, ]
  y = w * betaT + epsilon_alls[[1]]
  betaThat_ideals[nsim_ideal] = sum(y * w) / n
}
ideal_sq_err = mean((betaThat_ideals - betaT)^2)
log10(ideal_sq_err)