library(ggplot2)
library(tidyr)
library(dplyr)

sim_pois_coverage_single <- function(N_per_group, mu, b, sim_id) {
  group <- c(rep(0, N_per_group), rep(1, N_per_group))
  
  y <- rpois(2* N_per_group, lambda = exp(mu + b * group))
  data <- data.frame(y, group)
  
  cf_glm <- suppressMessages(confint(glm(y ~ group, data = data, family = "poisson"))[2, ])
  

  fit_rstanarm <- rstanarm::stan_glm(y ~ group, data = data, family = "poisson", cores = 1, chains = 2, refresh = 0)

  cf_rstanarm <- rstanarm::posterior_interval(fit_rstanarm, prob = 0.95)["group", ]

  data.frame(method = c("glm", "rstanarm"), sim_id, mu = mu, b = b, N_per_group = N_per_group,
             ci_low = unname(c(cf_glm[1], cf_rstanarm[1])), ci_high = unname(c(cf_glm[2], cf_rstanarm[2]))
  )
}



scenarios <- tidyr::crossing(N_per_group = c(1), mu = log(100), b = seq(0, 2, length.out = 10), sim_id = 1:100)
#scenarios <- tidyr::crossing(N_per_group = 4, mu = log(10), b = seq(0, 2, length.out = 2), sim_id = 1:10)

future::plan(future::multisession)

pois_coverage_df <- furrr::future_pmap_dfr(scenarios, sim_pois_coverage_single, .options = furrr::furrr_options(seed = TRUE, chunk_size = 40))

cache_dir <- "sim_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}
saveRDS(pois_coverage_df, file = file.path(cache_dir, "pois_coverage.rds"))


pois_coverage_df |> 
  filter(!is.na(ci_low), !is.na(ci_high)) |> 
  mutate(covered = b >= ci_low & b <= ci_high) |> 
  group_by(method, b, N_per_group) |> 
  summarise(coverage = mean(covered),
            n_covered = sum(covered),
            coverage_low = qbeta(0.025, n_covered, n() - n_covered + 1),
            coverage_high = qbeta(0.975, n_covered + 1 , n() - n_covered ),
            .groups = "drop"
  ) |> 
  ggplot()  + aes(x = b, ymin = coverage_low, y = coverage, ymax = coverage_high) +
  geom_hline(color = "orangered", yintercept = 0.95) +
  geom_ribbon(fill = "#888", alpha = 0.3) +
  geom_line() + facet_grid(method~N_per_group, labeller = "label_both") + 
  scale_x_continuous("True fold change") + theme(strip.text.y = element_text(size = 10))
