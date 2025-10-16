library(MASS)
library(rstanarm)
library(ggplot2)
theme_set(theme_minimal())

sim_nb_coverage_single <- function(N_per_group, mu, b, phi, sim_id = NULL) {
  group <- c(rep(0, N_per_group), rep(1, N_per_group))

  y <- rnbinom(2* N_per_group, mu = exp(mu + b * group), 
             size = phi)
  data <- data.frame(y, group)

  fit_glm.nb <- MASS::glm.nb(y ~ group, data = data)
  cf_glm.nb <- suppressMessages(confint(fit_glm.nb)["group", ])
  estimate_glm.nb <- coef(fit_glm.nb)["group"] 

  fit_rstanarm <- rstanarm::stan_glm.nb(y ~ group, data = data, refresh = 0, prior = NULL, prior_intercept = NULL)
  cf_rstanarm <- rstanarm::posterior_interval(fit_rstanarm, prob = 0.95)["group",]
  estimate_rstanarm <- coef(fit_rstanarm)["group"]
  
  res <- data.frame(method = c("glm.nb", "rstanarm"), true_b = b, 
                    estimate = c(estimate_glm.nb, estimate_rstanarm),
             ci_low = unname(c(cf_glm.nb[1], cf_rstanarm[1])), ci_high = unname(c(cf_glm.nb[2], cf_rstanarm[2])))
  if(!is.null(sim_id)) {
    res$sim_id <- sim_id
  }
  res
}
  
    
sim_nb_coverage_single(N_per_group = 4, mu = log(100), b = 1, phi = 0.5)

sim_nb_multiple <- function(N_sims, N_per_group, mu, b, phi) {
  res_list <- lapply(
    1:N_sims,
    FUN =
      \(sim_id) sim_nb_coverage_single(
        N_per_group = N_per_group,
        mu = mu,
        b = b,
        phi = phi,
        sim_id = sim_id
      )
  )
  do.call(rbind, res_list)
}

res_multiple <- sim_nb_multiple(5, N_per_group = 4, mu = log(100), b = 1, phi = 0.5)

ggplot(res_multiple) + aes(x = sim_id, y = estimate, ymin = ci_low, ymax = ci_high, color = method) +
  geom_hline(aes(yintercept = true_b), color = "blue") +
  geom_pointrange(position = position_dodge(width = 0.3))
