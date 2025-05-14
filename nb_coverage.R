library(ggplot2)
library(tidyr)
library(dplyr)

sim_nb_coverage_single <- function(N_per_group, mu, b, phi, sim_id, base_brm) {
  group <- c(rep(0, N_per_group), rep(1, N_per_group))
  
  y <- rnbinom(2* N_per_group, mu = exp(mu + b * group), 
               size = phi)
  data <- data.frame(y, group)
  
  cf_glm.nb <- suppressMessages(confint(MASS::glm.nb(y ~ group, data = data))[2, ])
  
  cf_gamlss <- tryCatch({
  fit_gamlss <- suppressMessages(gamlss::gamlss(y ~ group, family = "NBI"))
  confint(fit_gamlss)[2,]
  }, error = function(e) {return(c(NA,NA))})
  
  # Rarely we get problematic init which results in fitted phi -> infty and bad BFMI 
  # refitting fixes that (more informative prior on phi also would)
  for(i in 1:5) {
    fit_brm <- update(base_brm, newdata = data, cores = 1, chains = 2, future = FALSE, backend = "cmdstanr", refresh = 0)

    if(all(!is.na(rstan::get_bfmi(fit_brm$fit)))) {
      break
    }
  }
  bfmi_problem <- any(is.na(rstan::get_bfmi(fit_brm$fit)))
  cf_brm <- brms::fixef(fit_brm)["group", c(3,4)]
  # bfmi_problem <- FALSE
  # cf_brm <- c(NA,NA) 

  data.frame(method = c("glm.nb", "gamlss", "brms"), sim_id, mu = mu, b = b, phi = phi, N_per_group = N_per_group,
             ci_low = unname(c(cf_glm.nb[1], cf_gamlss[1], cf_brm[1])), ci_high = unname(c(cf_glm.nb[2], cf_gamlss[2], cf_brm[2])),
             bfmi_problem = c(FALSE, FALSE, bfmi_problem)
  )
}


# Construct base brms object to update
mu_0 <- log(100)
group_eff <- log(1.5)
true_phi <- 2.5
group <- c(rep(0, 4), rep(1, 4))
y <- rnbinom(8, mu = exp(mu_0 + group_eff * group), 
             size = true_phi)
prior <- c(
  brms::prior("", class = "Intercept")
)
base_brm <- brms::brm(y ~ group, data = data.frame(y, group), family = "negbinomial", prior = prior, backend = "cmdstanr")


scenarios <- tidyr::crossing(N_per_group = 4, mu = log(100), b = seq(0, 2, length.out = 10), phi = c(0.5, 1, 2.5, 10), sim_id = 1:1000)
#scenarios <- crossing(N_per_group = 4, mu = log(100), b = seq(0, 2, length.out = 5), phi = c(1, 2.5), sim_id = 1:2)

future::plan(future::multisession)

nb_coverage_df <- furrr::future_pmap_dfr(scenarios, \(...) sim_nb_coverage_single(..., base_brm = base_brm), .options = furrr::furrr_options(seed = TRUE, chunk_size = 40))

cache_dir <- "sim_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
saveRDS(nb_coverage_df, file = file.path(cache_dir, "nb_coverage.rds"))

