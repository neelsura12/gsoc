library(cmdstanr)
fp <- file.path("/Users/nshah/work/gsoc/models/lamw_normal.stan")
mod <- cmdstan_model(fp, force_recompile=T)

# tmp ---------------------------------------------------------------------

N <- 1000
mu <- 0
sigma <- 1
delta <- 1

x <- rnorm(N, mean=mu, sd=sigma)
x_hat=(x - mu)/sigma
y <- x_hat * exp(delta/2*x_hat^2)*sigma + mu;

mod_out <- mod$sample(
    data=list(N=N,
                y=yy,
                mu=mu,
                sigma=sigma),
    chains=2,
    init=5,
    adapt_delta=0.8,
    parallel_chains=2,
    iter_warmup=200,
    iter_sampling=200
)

mod_out$summary()

# ------------------

yy=LambertW::rLambertW(N, distname="normal",
                         theta=list(beta=c(mu, sigma), delta=delta))
