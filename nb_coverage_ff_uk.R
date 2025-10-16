library(ggplot2)
library(tidyr)
library(dplyr)
library(rstanarm)


sim_nb_coverage_single <- function(N_per_group, mu, b, phi, sim_id) {
  library(rstanarm)
  group <- c(rep(0, N_per_group), rep(1, N_per_group))
  
  y <- rnbinom(2* N_per_group, mu = exp(mu + b * group), 
               size = phi)
  data <- data.frame(y, group)
  
  cf_glm.nb <- suppressMessages(confint(MASS::glm.nb(y ~ group, data = data))[2, ])
  

  m <- glmmTMB::glmmTMB(y ~ group, family = glmmTMB::nbinom2, data = data)
  cf_glmmTMB <- confint(m, method = "profile", parm = "group", estimate = FALSE)
  
  fit_rstanarm <- rstanarm::stan_glm.nb(y ~ group, data = data, refresh = 0, prior = NULL, prior_intercept = NULL)
  cf_rstanarm <- rstanarm::posterior_interval(fit_rstanarm, prob = 0.95)["group",]

  data.frame(method = c("glm.nb", "glmmTMB", "rstanarm"), sim_id, mu = mu, b = b, phi = phi, N_per_group = N_per_group,
             ci_low = unname(c(cf_glm.nb[1], cf_glmmTMB[1], cf_rstanarm[1])), ci_high = unname(c(cf_glm.nb[2], cf_glmmTMB[2], cf_rstanarm[2]))
  )
}



scenarios <- tidyr::crossing(N_per_group = 4, mu = log(100), b = seq(0, 2, length.out = 10), phi = c(0.5, 1, 2.5, 10), sim_id = 1:100)
#scenarios <- crossing(N_per_group = 4, mu = log(100), b = seq(0, 2, length.out = 5), phi = c(1, 2.5), sim_id = 1:2)

future::plan(future::multisession)
#future::plan(future::sequential)

nb_coverage_df <- furrr::future_pmap_dfr(scenarios, sim_nb_coverage_single, .options = furrr::furrr_options(seed = TRUE, chunk_size = 40))

cache_dir <- "sim_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
saveRDS(nb_coverage_df, file = file.path(cache_dir, "nb_coverage_ff_uk.rds"))

