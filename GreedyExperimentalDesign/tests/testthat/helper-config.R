build_test_configs <- function() {
  n_values = c(16, 96, 192)
  p_values = c(1, 5, 10)
  B_values = c(2, 8)
  r_values = c(100, 200, 1000)
  seed_counter = 0
  configs = list()
  for (n in n_values) {
    for (B in B_values) {
      if (n %% B != 0) {
        next
      }
      n_B = n / B
      if (n_B < 4) {
        next
      }
      prop_T_values = unique(c(1 / n_B, floor(n_B / 2) / n_B))
      prop_T_values = prop_T_values[prop_T_values > 0 & prop_T_values < 1]
      for (p in p_values) {
        for (r in r_values) {
          for (prop_T in prop_T_values) {
            seed_counter = seed_counter + 1
            configs[[length(configs) + 1]] = list(
              n = n,
              p = p,
              r = r,
              B = B,
              prop_T = prop_T,
              seed = seed_counter
            )
          }
        }
      }
    }
  }
  configs
}
