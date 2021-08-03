library(cmdstanr)
library(LambertW)
library(posterior)

N <- 1000
delta <- 0.1

# h normal
s <- rpois(N, lambda=3)
sigma <- 1
x <- 3*s + 0.5 + rnorm(N, sd=sigma)

mu <- 3*s + 0.5
x_hat = (x - mu)/sigma
y <- x_hat * exp(delta/2*x_hat^2)*sigma + mu;

# classic
# y = norm(beta*x + alpha, error_sd)

# what we're doing
# y=LambertW x Fx X 
 # X = norm(beta*s+alpha, error_sd)
 # Y = u * exp(delta/2 * u^2) where u = (x-mu)/sigma

fp <- file.path("/Users/nshah/work/gsoc/models/lambertw_normal_h_transform.stan")
mod <- cmdstan_model(fp, force_recompile = F)

# works for all distributions
mod_out <- mod$sample(data=list(N=N, y=y, s=s), parallel_chains=4)
mod_out$summary()

quantile(as_draws_df(mod_out$draws())$new_y, c(0.05,0.25,0.5,0.75,0.95))
quantile(y, c(0.05,0.25,0.5,0.75,0.95))