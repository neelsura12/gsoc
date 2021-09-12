library(cmdstanr)
library(LambertW)


N <- 100
mu <- 0
sigma <- 1

y <- LambertW::rLambertW(N, "normal", theta=list(beta=c(mu, sigma), 
                                            gamma=0.1,
                                            alpha=1,
                                            delta=c(0,0)))

# fit stan program ------------------------------------------------------------

fp <- file.path(paste(getwd(), "/week36/lambertw_normal_s_transform.stan", sep=""))
mod <- cmdstan_model(fp, force_recompile = F)

mod_out <- mod$sample(data=list(N=N, y=y, mu=mu, sigma=sigma), parallel_chains=4)
mod_out$summary()