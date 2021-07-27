library(cmdstanr)
library(posterior)
library(invgamma)
library(LambertW)

c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

fp <- file.path("/Users/nshah/Documents/GitHub/scratch/models/lamw_normal_explicit.stan")
fit_model <- cmdstan_model(fp, force_recompile=T)

# -----------------------------------
N <- 100  # num observations
R <- 1000 # num draws from joint
S <- 500  # num SBC iters

prior_mu <- rnorm(S, 0, 1)
prior_sigma <- abs(rnorm(S, sqrt(pi/2)))
prior_delta <- rexp(S, rate=1)

sbc_rank_delta=c()
sbc_rank_mu=c()
sbc_rank_sigma=c()
idx = seq(1, 4000 - 8, 8)
for (i in 1:S) {
  y_simu <- LambertW::rLambertW(N, distname="normal", 
                                theta=list(beta=c(prior_mu[i], prior_sigma[i]), 
                                           delta=prior_delta[i], 
                                           alpha=1))
  mod_out <- fit_model$sample(data=list(N=N,
                                        y=y_simu),
                              iter_sampling=R,
                              parallel_chains=4,
                              show_messages=F,
                              adapt_delta=0.8,
                              iter_warmup=R,
                              refresh=0)
  posterior_samples <- as_draws_df(mod_out$draws())
  
  sbc_rank_delta=c(sbc_rank_delta, sum(prior_delta[i] < posterior_samples$delta[idx]))
  sbc_rank_mu=c(sbc_rank_mu, sum(prior_mu[i] < posterior_samples$mu[idx]))
  sbc_rank_sigma=c(sbc_rank_sigma, sum(prior_sigma[i] < posterior_samples$sigma[idx]))
  if (i %% 100 == 0) {
  	cat("Res:", i, "\n")
  }
}

# ------------------------------------------------------------------------------
sbc_hist <- hist(sbc_rank_mu/S, seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC mu", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)

sbc_hist <- hist(sbc_rank_sigma/S, seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC sigma", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)

sbc_hist <- hist(sbc_rank_delta/S, seq(0, 1, 0.05), col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="SBC delta", xlab="Prior Rank", yaxt='n', ylab="", col=c_dark, border=c_dark_highlight)
