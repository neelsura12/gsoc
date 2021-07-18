library(cmdstanr)
set_cmdstan_path("/Users/pinkney/.cmdstanr/cmdstan-2.27.0")
library(posterior)
library(invgamma)
library(LambertW)

c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

fp <- file.path("~../lamw_normal_joint.stan")
fit_model <- cmdstan_model(fp, force_recompile=F)

# -----------------------------------
N <- 100  # num observations
R <- 1000 # num draws from joint
S <- 100  # num SBC iters

prior_mu <- rnorm(N, 0, 1)
prior_sigma <- rinvgamma(N, shape=2, scale=.1)
prior_delta <- rexp(N, rate=1)

sbc_rank_delta=c()
sbc_rank_mu=c()
sbc_rank_sigma=c()
for (i in 1:S) {
  y_simu <- LambertW::rLambertW(N, distname="normal", 
                                theta=list(beta=c(prior_mu[i], prior_sigma[i]), 
                                           delta=prior_delta[i], 
                                           alpha=1))
  mod_out <- fit_model$sample(data=list(N=N,
                                        y=y_simu),
                              iter_sampling=R,
                              parallel_chains = 4,
                              show_messages = F,
                              adapt_delta = 0.8,
                              iter_warmup=R,
                              refresh=0)
  posterior_samples <- as_draws_df(mod_out$draws())
  
  sbc_rank_delta=c(sbc_rank_delta, sum(prior_delta[i] < posterior_samples$delta))
  sbc_rank_mu=c(sbc_rank_mu, sum(prior_mu[i] < posterior_samples$mu))
  sbc_rank_sigma=c(sbc_rank_sigma, sum(prior_sigma[i] < posterior_samples$sigma))
}

# ------------------------------------------------------------------------------
sbc_hist <- hist(sbc_rank_mu/(4*R), seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC mu", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)

sbc_hist <- hist(sbc_rank_sigma/(4*R), seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC sigma", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)

sbc_hist <- hist(sbc_rank_delta/(4*R), seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC delta", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)
