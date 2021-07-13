library(cmdstanr)
set_cmdstan_path("/Users/nshah/cmdstan")
fp <- file.path("/Users/nshah/Documents/GitHub/scratch/models/lamw_tdist.stan")

# ------------------

N <- 1000
mu <- 0
nu <- 10
sigma <- sqrt(nu/(nu - 2))
delta <- 1/3

x <- rt(N, df=nu) + mu
x_hat <- (x - mu)/sigma
y <- x_hat * exp(delta/2 * x_hat^2) * sigma + mu;

mod <- cmdstan_model(fp, force_recompile = T)
mod_out <- mod$sample(
    data = list(N = N, y = y, nu = nu, mu = mu, sigma = sigma),
    chains = 2,
    init = 5,
    adapt_delta = 0.8,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 200
)

mod_out$summary()

#  ------------------

yy <- LambertW::rLambertW(N, distname="t", theta=list(beta=c(location=mu, scale=sigma, df=nu), delta=delta))
LambertW::MLE_LambertW(yy, distname="t", type="h", theta.fixed=list(beta=c(location=mu, scale=sigma, df=nu), alpha=1))

