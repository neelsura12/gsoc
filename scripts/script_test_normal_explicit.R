library(cmdstanr)
fp <- file.path("/Users/nshah/Documents/GitHub/scratch/models/lamw_normal_explicit.stan")
mod <- cmdstan_model(fp, force_recompile=T)

# tmp ---------------------------------------------------------------------

N <- 1000
mu <- 0
sigma <- 1
delta <- 1

y = LambertW::rLambertW(N, theta=list(beta=c(mu, sigma), delta=delta, alpha=1), distname="normal")

mod_out <- mod$sample(
    data=list(N=N, y=y),
    parallel_chains=4
)

mod_out$summary()