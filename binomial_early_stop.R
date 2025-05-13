library(dplyr)
library(tidyr)
sim_bayesian_early_stop_single <- function(true_prop, low_stop, high_stop, width_stop, max_steps) {
  step_size <- 100
  steps <- numeric(max_steps)
  start <- 1
  N_success_prev <- 0
  while(start <= max_steps) {
    to_sim <- min(step_size, max_steps - start + 1)
    y <- rbinom(to_sim, size = 1, prob = true_prop)
    N <- (start - 1) + 1:to_sim
    N_success <- N_success_prev + cumsum(y)
    
    ci_low <- qbeta(0.025, N_success + 1, N - N_success + 1)
    ci_high <- qbeta(0.975, N_success + 1, N - N_success + 1)
    stop <- FALSE
    if(!is.na(low_stop)) {
      stop <- stop | ci_high < low_stop
    }
    if(!is.na(high_stop)) {
      stop <- stop | ci_low > high_stop
    }
    if(!is.na(width_stop)) {
      stop <- stop | ci_high - ci_low < width_stop
    }
    if(any(stop)) {
      stop_index <- which(stop)[1]
      return(
        data.frame(true_prop, low_stop, high_stop, width_stop, stopped = TRUE, N = N[stop_index], N_success = N_success[stop_index], ci_low = ci_low[stop_index], ci_high = ci_high[stop_index])
      )
    } 
    start <- start + to_sim
    N_success_prev <- N_success_prev + sum(y) 
  }
  data.frame(true_prop, low_stop, high_stop, width_stop, stopped = FALSE, N = max_steps, N_success = N_success_prev, ci_low = ci_low[to_sim], ci_high = ci_high[to_sim])
}


sim_bayesian_early_stop <- function(true_props, N_sims, low_stop, high_stop, width_stop, max_steps) {
  cache_dir <- "sim_cache"
  if(!dir.exists(cache_dir)) {
    dir.create(cache_dir)
  }
  
  cache_hash <- rlang::hash(list(
    true_props, N_sims, low_stop, high_stop, width_stop, max_steps
  ))
  cache_file <- file.path(cache_dir, paste0("early_stop_", cache_hash,".rds"))
  
  if(file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  future::plan(future.mirai::mirai_multisession)
  
  res <- furrr::future_map_dfr(
    true_props,
    \(true_prop) {
      purrr::map_dfr(1:N_sims, \(x) { sim_bayesian_early_stop_single(true_prop = true_prop, low_stop = low_stop, high_stop = high_stop, width_stop = width_stop, max_steps = max_steps)})  
    },
    .options = furrr::furrr_options(seed = TRUE, chunk_size = 1)
  )
  saveRDS(res, file = cache_file)
  
  res
}

plot_coverage_stop <- function(coverage_stop) {
  coverage_stop |> 
    mutate(covered = true_prop >= ci_low & true_prop <= ci_high) |> 
    group_by(true_prop) |> 
    summarise(coverage = mean(covered),
              n_covered = sum(covered),
              coverage_low = qbeta(0.025, n_covered, n() - n_covered + 1),
              coverage_high = qbeta(0.975, n_covered + 1 , n() - n_covered )
              ) |> 
    ggplot()  + aes(x = true_prop, ymin = coverage_low, y = coverage, ymax = coverage_high) +
      geom_hline(color = "orangered", yintercept = 0.95) +
      geom_ribbon(fill = "#ccc") +
      geom_line() +
      scale_x_continuous("True proportion")
}



coverage_stop_05 <- sim_bayesian_early_stop(true_props = seq(0.05, 0.95, by = 0.01), N_sims = 1000, low_stop = 0.5, high_stop = 0.5, width_stop = NA, max_steps = 10000)


coverage_stop_04_06 <- sim_bayesian_early_stop(true_props = seq(0.05, 0.95, by = 0.01), N_sims = 1000, low_stop = 0.4, high_stop = 0.6, width_stop = NA, max_steps = 10000)

coverage_stop_width_01 <- sim_bayesian_early_stop(true_props = seq(0.05, 0.95, by = 0.01), N_sims = 1000, low_stop = NA, high_stop = NA, width_stop = 0.1, max_steps = 10000)

