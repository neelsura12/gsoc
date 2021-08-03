library(cmdstanr)
library(LambertW)

N <- 1000
delta_left <- 0.1
delta_right <- 0.1

# hh normal
mu <- 0
sigma <- 1

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/models/lambertw_normal_hh_transform.stan")
mod <- cmdstan_model(fp, force_recompile = F)
y <- LambertW::rLambertW(
    N, distname="normal", theta=list(beta=c(mu, sigma), 
                                     delta=c(delta_left, delta_right)))

# works for all disributions
mod_out <- mod$sample(data=list(N=N, y=y), parallel_chains=4)
mod_out$summary()
