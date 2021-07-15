library(cmdstanr)
fp <- file.path("/Users/nshah/work/gsoc/models/lamw_normal.stan")
fit_model <- cmdstan_model(fp, force_recompile = T)

# ------------------

N <- 100
prior_mu <- rnorm(N, 0, 1)
prior_sd <- sqrt(1/rgamma(N, shape=1))
prior_delta <- rexp(N, rate=1)

for (i in 1:N) {
    y_tmp <- LambertW::rLambertW(N, distname = "normal", 
                                 theta=list(beta=c(prior_mu[i], prior_sd[i]), 
                                            delta = prior_delta[i]))
}

# tmp ---------------------------------------------------------------------

tryCatch({
    registerDoParallel(makeCluster(detectCores()))
    
    simu_list <- t(data.matrix(data.frame(simu_lambdas, simu_ys)))
    
    ensemble_output <- foreach(simu=simu_list,
                               .combine='cbind') %dopar% {
                                   simu_lambda <- simu[1]
                                   simu_y <- simu[2:(N + 1)];
                                   
                                   # Fit the simulated observation
                                   input_data <- list("N" = N, "y" = simu_y)
                                   
                                   capture.output(library(rstan))
                                   capture.output(fit <- sampling(fit_model, data=input_data, seed=4938483))
                                   
                                   # Compute rank of prior draw with respect to thinned posterior draws
                                   sbc_rank <- sum(simu_lambda < extract(fit)$lambda[seq(1, N - 2, 2)])
                                   
                                   c(sbc_rank)
                               }
}, finally={ stopImplicitCluster() })

# Check SBC histogram
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="Lambda", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)
