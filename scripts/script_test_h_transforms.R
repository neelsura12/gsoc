library(cmdstanr)
library(LambertW)

N <- 1000
delta <- 1

# h normal
mu <- 0
sigma <- 1

fp <- file.path("/Users/nshah/work/gsoc/models/lambertw_normal_h_transform.stan")
mod <- cmdstan_model(fp, force_recompile = F)
y <- LambertW::rLambertW(
    N, distname="normal", theta=list(beta=c(mu, sigma), delta=delta))

# h exponential
lambda <- 1

fp <- file.path("/Users/nshah/work/gsoc/models/lambertw_exponential_h_transform.stan")
mod <- cmdstan_model(fp, force_recompile = F)
y <- LambertW::rLambertW(
    N, distname="exp", theta=list(beta=c(rate=lambda), delta=delta))

# works for all disributions
mod_out <- mod$sample(data=list(N=N, y=y), parallel_chains=4)
mod_out$summary()
