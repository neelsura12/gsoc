library(cmdstanr)
library(posterior)

# h
fp <- file.path("~/work/gsoc/models/lambertw_normal_h_rng.stan")
mod <- cmdstan_model(fp, force_recompile=F)
mod_out <- mod$sample(data=list(N=100, mu=0, sigma=1, delta=0.1), fixed_param=T)

# hh
fp <- file.path("~/work/gsoc/models/lambertw_normal_hh_rng.stan")
mod <- cmdstan_model(fp, force_recompile=F)
mod_out <- mod$sample(
    data=list(N=100, mu=0, sigma=1, delta_left=0.1, delta_right=0), 
    fixed_param=T)

# s
fp <- file.path("~/work/gsoc/models/lambertw_normal_s_rng.stan")
mod <- cmdstan_model(fp, force_recompile=F)
mod_out <- mod$sample(
    data=list(N=100, mu=0, sigma=1, gamma=-1), 
    fixed_param=T)

# plot (works for all types)
hist(as_draws_df(mod_out$draws())$x)
hist(as_draws_df(mod_out$draws())$y)
