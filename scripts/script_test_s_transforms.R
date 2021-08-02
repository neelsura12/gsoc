library(cmdstanr)
library(LambertW)

N <- 1000
gamma <- 0.5

# s exponential
lambda <- 1

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/models/lambertw_exponential_s_transform.stan")
mod <- cmdstan_model(fp, force_recompile = F)
y <- LambertW::rLambertW(
    N, distname="exp", theta=list(beta=c(lambda), gamma=gamma))

# works for all disributions
mod_out <- mod$sample(data=list(N=N, y=y), parallel_chains=4)
mod_out$summary()
