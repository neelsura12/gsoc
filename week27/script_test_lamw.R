library(cmdstanr)
set_cmdstan_path('/Users/nshah/cmdstan')
fp <- file.path("./lambertw.stan")

# ------------------

N <- 1000
mu <- 0
sigma <- 1
delta <- 1

x <- rnorm(N, mean = mu, sd = sigma)
x_hat <- (x - mu)/sigma
y <- x_hat * exp(delta/2 * x_hat^2) * sigma + mu;

mod <- cmdstan_model(fp, force_recompile = T)
mod_out <- mod$sample(
    data = list(N = N, y = y),
    chains = 2,
    init = 5,
    adapt_delta = 0.8,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 200
)

mod_out$summary()

# ------------------

library(LambertW)
yy = LambertW::rLambertW(N, beta = c(mu, sigma), delta = delta, distname = "normal")
