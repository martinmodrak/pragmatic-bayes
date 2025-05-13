library(ggplot2)
clopper_pearson_ci_fun <- function(N, N_success, alpha = 0.05) {
  list(
    low = qbeta(alpha / 2, N_success, N - N_success + 1),
    high = qbeta(1 - alpha / 2, N_success + 1, N - N_success)
  )
}

bayesian_uniform_ci_fun <- function(N, N_success, alpha = 0.05) {
  list(
    low = qbeta(alpha / 2, N_success + 1, N - N_success + 1),
    high = qbeta(1 - alpha / 2, N_success + 1, N - N_success  + 1)
  )
}

bayesian_uniform_hpdi_fun <- function(N, N_success, alpha = 0.05) {
  low <- numeric(length(N_success))
  high <- numeric(length(N_success))
  for(i in seq_along(N_success)) {
    if(N_success[i] == 0) {
      low[i] <- 0
      high[i] <- qbeta(1 - alpha, N_success[i] + 1, N - N_success[i] + 1)
    } else if(N_success[i] == N) {
      low[i] <- qbeta(alpha, N_success[i] + 1, N - N_success[i] + 1)
      high[i] <- 1
    } else {
      hdi_ <- HDInterval::hdi(qbeta, shape1 = N_success[i] + 1, shape2 = N - N_success[i] + 1)
      low[i] <- hdi_[1]
      high[i] <- hdi_[2]
    }
  }
  list(low = low, high = high)
}


coverage_single <- function(N, proportion, ci, alpha = 0.05) {
  weights <- dbinom(0:N, size = N, prob = proportion)
  coverage_all <- ci$low <= proportion & ci$high >= proportion
  
  sum(weights * coverage_all)
}

coverage_df <- function(N, props_to_test, ci_fun, label, alpha = 0.05) {
  ci <- ci_fun(N, 0:N, alpha = alpha)
  
  data.frame(
    N = N,
    props_to_test = props_to_test,
    label = label,
    coverage = purrr::map_dbl(props_to_test, .f = \(x) coverage_single(N = N_to_test, proportion = x, ci = ci))
  )
}

props_to_test <- seq(0, 1, by = 0.001)
N_to_test <- 20

coverage_df_bayesian_uniform <- 
  coverage_df(N_to_test, props_to_test, label = "Bayesian CI\n(uniform prior)", ci_fun = bayesian_uniform_ci_fun)

coverage_df_bayesian_hdpi <- 
  coverage_df(N_to_test, props_to_test, label = "Bayesian HPDI\n(uniform prior)", ci_fun = bayesian_uniform_hpdi_fun)

coverage_df_cp <- 
  coverage_df(N_to_test, props_to_test, label = "Clopper-Pearson", ci_fun = clopper_pearson_ci_fun)

coverage_for_plot <- rbind(
  coverage_df_bayesian_uniform,
  coverage_df_bayesian_hdpi,
  coverage_df_cp
)

coverage_plot <- function(coverage_for_plot) {
  coverage_for_plot |> ggplot() + aes(x = props_to_test, y = coverage) +
    geom_hline(color = "orangered", yintercept = 0.95) + 
    geom_line() + 
    facet_wrap(~label) +
    coord_cartesian(ylim = c(0.9, 1))
  
}

print(
  coverage_plot(coverage_for_plot)
)
