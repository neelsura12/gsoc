library(cmdstanr)
fp <- file.path("/Users/nshah/work/gsoc/models/lamw_exponential.stan")
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
yy <- LambertW::rLambertW(N, distname="exp", theta=list(beta=c(rate=lambda), delta=delta))

mod_out <- mod$sample(
    data = list(N = N, y = y, lambda=lambda),
    chains = 2,
    init = 5,
    adapt_delta = 0.8,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 200
)

mod_out$summary()

# ------------------

yy <- LambertW::rLambertW(N, distname="exp", theta=list(beta=c(rate=lambda), delta=delta))
LambertW::MLE_LambertW(yy, distname="exp", type="h", theta.fixed=list(beta=c(rate=lambda), alpha=1))

Neg binomial ~ LambertW x Exp (lambda, delta)

# ------------------


# figure this for weibull
ytrue <- rweibull(100, scale=1, shape=0.5) 
mu <- mean(ytrue)

# figure out delta
mle <- LambertW::MLE_LambertW(ytrue, distname="exp", type="h", theta.fixed=list(beta=c(rate=mu), alpha=1))

# simulate lambertW x F_X where F_X is exp
ylamw <- LambertW::rLambertW(100, distname="exp", theta=list(beta=c(rate=mu), delta=mle$tau["delta"]))

# ---
quantile(ytrue, c(0.025,0.25,0.5,0.75,0.975))
quantile(ylamw, c(0.025,0.25,0.5,0.75,0.975))
