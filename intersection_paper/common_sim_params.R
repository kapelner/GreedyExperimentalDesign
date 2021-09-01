#https://console.aws.amazon.com/cost-management/home?#/custom?reportType=CostUsage&chartStyle=Stack&groupBy=InstanceType&forecastTimeRangeOption=None&hasBlended=false&hasAmortized=false&excludeDiscounts=true&usageAs=usageQuantity&excludeCategorizedResources=false&excludeTaggedResources=false&excludeForecast=false&reportName=&timeRangeOption=Last6Months&startDate=2020-12-01&endDate=2021-05-31&granularity=Monthly&filter=%5B%5D&isTemplate=true




Xalls = list()
for (covariate_distribution in covariate_distributions){
  Xalls_cov = list()
  for (nx in 1 : nX){
    if (covariate_distribution == "uniform"){
      X = matrix(runif(n * max(ps), -sqrt(3), sqrt(3)), nrow = n)
    } else if (covariate_distribution == "normal"){
      X = matrix(rnorm(n * max(ps)), nrow = n)
    } else if (covariate_distribution == "exponential"){
      X = matrix(rexp(n * max(ps)), nrow = n)
    }
    Xalls_cov[[nx]] = X
    
    # Xtemp = matrix(rnorm(n * p, 0, 1), nrow = n)
    # Xtemp = data.matrix(MASS::Pima.tr[1 : n, 1 : p])
    # Xalls[[nx]] = apply(Xtemp, 2, function(xj){(xj - mean(xj)) / sd(xj)}) #stdize
    # rm(Xtemp)
  }
  Xalls[[covariate_distribution]] = Xalls_cov
}

epsilon_alls = list()
for (neps in 1 : nEPS){
  epsilon_alls[[neps]] = rnorm(n, 0, sigma_e)
}

all_betas_and_correlations = list()
all_betas_and_correlations$zero_effects = list()
all_betas_and_correlations$linear_uniform = list(
  beta_x1 = 3,
  beta_x2 = 3
)
all_betas_and_correlations$one_quadratic = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1
)
all_betas_and_correlations$quadratics_and_interaction = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1
)
all_betas_and_correlations$purely_nonlinear = list(
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1
)
all_betas_and_correlations$linear_uniform_plus_hetero_tx_effect = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1_T = 1,
  beta_x2_T = 1
)
all_betas_and_correlations$quadratics_and_interaction_plus_hetero_tx_effect = list(
  beta_x1 = 3,
  beta_x2 = 3,
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 1,
  beta_x2_T = 1
)
all_betas_and_correlations$purely_nonlinear_plus_hetero_tx_effect = list(
  beta_x1sq = 1,
  beta_x2sq = 1,
  beta_x1_x2 = 1,
  beta_x1_T = 1,
  beta_x2_T = 1
)
all_betas_and_correlations$purely_nonlinear_sin = list(
  beta_x1_sin = 10, #to make lots of wiggles
  beta_x2_sin = 10, #to make lots of wiggles
  beta_x1_x_2_sin = 5 #to make lots of wiggles
)
all_betas_and_correlations$crazy_exponential = list(
  beta_x1_exp = 1,
  beta_x2_exp = 1,
  beta_x1_x_2_exp = 1
)

