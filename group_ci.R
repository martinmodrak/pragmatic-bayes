# A relatively stupid exploration based on https://stats.stackexchange.com/a/662565/73129
# Need to fix/handle convergence errors

library(MASS)
library(dplyr)
library(ggplot2)
library(brms)
options(brms.backend = "cmdstanr", mc.cores = 4)

mu_0 <- log(100)
group_eff <- log(1.5)
true_theta <- 2.5
group <- c(rep(0, 4), rep(1, 4))
# proof of concept
y <- rnegbin(8, mu = exp(mu_0 + group_eff * group), 
             theta = true_theta)
confint(glm.nb(y ~ group))[2, ]

prior <- c(
  prior("", class = "Intercept")
)
fit_brm <- brms::brm(y ~ group, data = data.frame(y, group), family = "negbinomial", prior = prior)
fixef(fit_brm)["group", c(3,4)]
# simulation
nsim <- 100
# calculates relatively fast, but spams messages
res <- t(replicate(nsim, {
  y <- rnegbin(8, mu = exp(mu_0 + group_eff * group), 
               theta = true_theta)
  cf_freq <- suppressMessages(confint(glm.nb(y ~ group))[2, ])
  for(i in 1:3) {
    fit_brm <- brms::brm(y ~ group, data = data.frame(y, group), family = "negbinomial", 
                         refresh = 0, 
                         init = 0.1,
                         prior = prior)
    
    if(all(!is.na(rstan::get_bfmi(fit_brm$fit)))) {
      break
    }
  }
  if(any(is.na(rstan::get_bfmi(fit_brm$fit)))) {
    stop("BFMI")
  }
  cf_brm <- fixef(fit_brm)["group", c(3,4)]
  c(cf_freq, cf_brm)
})
)

dat <- data.frame(res, group_eff, i = 1:nsim)
# just for fun
ggplot(dat, aes(y = rank(`X2.5..`), x = group_eff)) +
  geom_linerange(aes(xmin = `X2.5..`, xmax = `X97.5..`, 
                     color = group_eff < X97.5.. & group_eff > X2.5..)) +
  geom_point(color = "red")  

ggplot(dat, aes(y = rank(Q2.5), x = group_eff)) +
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5, 
                     color = group_eff <  Q97.5 & group_eff > Q2.5)) +
  geom_point(color = "red")  

# about 86% coverage for the 95%-CI
dat %>% filter(group_eff < X97.5.. & group_eff > X2.5..) %>%  count()/nsim
dat %>% filter(group_eff < Q97.5 & group_eff > Q2.5) %>%  count()/nsim
