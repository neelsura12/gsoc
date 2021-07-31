library(cmdstanr)
fp <- file.path("/Users/nshah/work/gsoc/models/lamw_exponential_joint.stan")
mod <- cmdstan_model(fp, force_recompile = T)

# ------------------

N <- 1000
lambda <- 1
mu <- 1/lambda
sigma <- 1/lambda
delta <- 5

x <- rexp(N, rate=lambda)
x_hat <- (x - mu)/sigma
y <- x_hat * exp(delta/2 * x_hat^2) * sigma + mu;

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

yy <- LambertW::rLambertW(N, distname="t", theta=list(beta=c(location=mu, scale=sigma, df=nu), delta=delta))
LambertW::MLE_LambertW(yy, distname="t", type="h", theta.fixed=list(beta=c(location=mu, scale=sigma, df=nu), alpha=1))

