library(cmdstanr)

N <- 1000
lambda <- 3
b <- 3
a <- 0.5
x <- rpois(N, lambda=lambda)
mu <- b*x + a
sd <- 1.5
nu <- 2*sd^2/(sd^2 - 1)
delta <- 1/3

# regression with normal error

y <- mu + rnorm(N, mu, sd)

qqnorm(y - mu, pch = 1, frame = FALSE, main="QQ plot: norm errors")
qqline(y - mu, col = "steelblue", lwd = 2)

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/week31/regression_w_normal_error.stan")
mod <- cmdstan_model(fp, force_recompile = F)

mod_out <- mod$sample(data=list(N=N, y=y, x=x), parallel_chains=4)
mod_out$summary()

# regression with t error

y <- mu + rt(N, nu)

qqnorm(y - mu, pch=1, frame=FALSE, ylim=c(-6,6), xlim=c(-3,3), main="QQ plot: t-dist errors")
qqline(y - mu, col="steelblue", lwd=2)

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/week31/regression_w_t_error.stan")
mod <- cmdstan_model(fp, force_recompile = F)

mod_out <- mod$sample(data=list(N=N, y=y, x=x), parallel_chains=4)
mod_out$summary()

# regression with non-linear transform

x <- mu + rnorm(N)
x_hat <- (x - mu) / sigma
y <-x_hat * exp(delta/2*x_hat^2)*sigma + mu

qqnorm(y - mu, pch=1, frame=FALSE, ylim=c(-6,6), xlim=c(-3,3), main="QQ plot: non-linear transform")
qqline(y - mu, col="steelblue", lwd=2)

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/week31/regression_w_non_linear_error.stan")
mod <- cmdstan_model(fp, force_recompile = F)

mod_out <- mod$sample(data=list(N=N, y=y, x=x), parallel_chains=4)
mod_out$summary()

# works on all 

quantile(as_draws_df(mod_out$draws())$new_y, c(0.05,0.25,0.5,0.75,0.95))
quantile(y, c(0.05,0.25,0.5,0.75,0.95))