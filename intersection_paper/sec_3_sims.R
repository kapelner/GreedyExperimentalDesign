####population model sims
nX = 50
nEPS = 50
nR = 500

source("common_sim_params.R")

exp_settings = data.frame(matrix(NA, nrow = nx * nEPS * length(ps), ncol = 3))
names(exp_settings) = c("p", "nx", "neps")
exp_settings$p =    rep(ps, each = nX * nEPS)
exp_settings$nx =   rep(rep(1 : nX, each = nEPS), length(ps))
exp_settings$neps = rep(1 : nEPS, nX * length(ps))
exp_settings

filename = "sec_S3"

source("common_sims.R")
