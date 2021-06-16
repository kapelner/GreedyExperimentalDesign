
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




pacman::p_load(stringr, dplyr)
res = tibble(design = factor(), model = factor(), p = numeric(), nx = integer(), neps = integer(), log10_squared_error = numeric())
for (filename in dir()){
  if (!str_detect(filename, "^sec_S3.*\\.RData$")){
    next
  }
  load(filename)
  inner_res = inner_res %>% group_by(design, model, p, nx, neps) %>% summarize(log10_squared_error = mean(log10(squared_error)))
  res = bind_rows(res, inner_res)
  # all_ltgr_data_by_file[[csv_file]] = read_csv(paste0("data/", csv_file))
  # csv_file_df = cbind(TRUE, str_replace(csv_file, ".csv", ""), read_csv(paste0("data/", csv_file)))
  print(filename)
  print(memory.size())
  Sys.sleep(0.05)
}



