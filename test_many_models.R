library(cmdstanr)
library(LambertW)

# normal h transform -----------------------------------------------------------

N <- 100
mu <- 0
sigma <- 1

y <- LambertW::rLambertW(N, "normal", theta=list(beta=c(mu, sigma), gamma=0.1, alpha=1, delta=c(0,0)))

fp <- file.path(paste(getwd(), "/week37/lambertw_normal_h.stan", sep=""))
mod <- cmdstan_model(fp, force_recompile = F)

mod_out <- mod$sample(data=list(N=N, y=y, mu=mu, sigma=sigma), parallel_chains=4)
mod_out$summary()
