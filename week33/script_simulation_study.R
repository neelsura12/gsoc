library(cmdstanr)
library(posterior)
library(LambertW)
library(xtable)

# setup -----------------------------------------------------------------------

epsilon <- function(N, sigma, delta_left, delta_right) {
    u <- rnorm(N)
    ifelse(u <= 0, 
           u*exp(delta_left/2*u^2)*sigma, 
           u*exp(delta_right/2*u^2)*sigma)
}

N <- 1000; mu_x <- 1; alpha <- 1; beta <- 3; sigma <- 3/2

x <- rnorm(N, mu_x, 1)

y1 <- alpha + beta*x + epsilon(N, sigma, 0, 1/3)
y2 <- alpha + beta*x + epsilon(N, sigma, 2/3, 1/3)
y <- c(y1, y2)
dim(y) <- c(N, 2)

# fit stan program ------------------------------------------------------------

fp <- file.path(paste(getwd(), "/week33/regression_lambertw_normal_hh.stan", sep=""))
mod <- cmdstan_model(fp, force_recompile = F)

for (i in 1:1) {
    mod_out <- mod$sample(data=list(N=N, y=y[,i], x=x), parallel_chains=4)
    print(xtable(mod_out$summary()[2:6,], type = "latex"), file=paste(getwd(), "/y", i, ".tex", sep=""))
}

# Plot distribution of y ------------------------------------------------------

# error histograms
par(mfrow = c(2, 2),
    cex = 0.6,
    mar = c(0, 0, 0, 0), 
    oma = c(3, 2, 0.5, 0.5),
    tcl = -0.25,
    mgp = c(2, 0.6, 0))
for (i in 1:2) {
    hist(y[,i], axes = T, freq = F, main=c(), yaxt="n")
    padj <- ifelse(i %in% c(3, 4), 1.5, 1)
    mtext(paste("Y", i, sep=""), side = 3, line = -1, adj = 0.025, padj = padj, cex = 2, col = "grey40")
    box(col = "grey60")
}
# error qqplots
for (i in 1:2) {
    qqnorm(y[,i], pch = 1, frame = FALSE, main=c())
    qqline(y[,i], col = "steelblue", lwd = 2)
    mtext(paste("Y", i, sep=""), side = 3, line = -1, adj = 0.025, padj = 1.5, cex = 2, col = "grey40")
    box(col = "grey60")
}

# Posterior retrodictive checks -----------------------------------------------

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

y_new <- as_draws_df(mod_out$draws('y_new'))$y_new

min_y = as.integer(min(y1, y_new) - 1)
max_y = as.integer(max(y1, y_new) + 2)
obs_counts <- hist(y1, breaks=(min_y:max_y)-0.5, plot=FALSE)$counts
B <- length(obs_counts) - 1

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)
x <- x - min_y
pad_obs <- do.call(cbind, lapply(idx, function(n) obs_counts[n + 1]))

counts <- sapply(1:4000, function(n) hist(y_new, breaks=(min_y:max_y)-0.5, plot=FALSE)$counts)
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:(B + 1), function(b) quantile(counts[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n + 1]))

plot(1, type="n", main="Posterior Retrodictive Check", xlim=c(-6, 32), xlab="y", ylim=c(0, max(c(obs_counts, cred[9,]))), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])), col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])), col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])), col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])), col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

lines(x, pad_obs, col="white", lty=1, lw=2.5)
lines(x, pad_obs, col="black", lty=1, lw=2)

# Bijective transform  --------------------------------------------------------

par(mfrow = c(1, 2))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
par(pch = 20)

delta <- 1/3

lamw_transform_h <- function(x, delta) {
    y <- x * exp(delta/2 * x^2)
}
lamw_untransform_h <- function(y, delta) {
    delta <- max(delta, 1e-16)
    x <- sign(y) * sqrt(W(delta * y^2) / delta)
}

lamw_transform_hh <- function(x, delta_left, delta_right) {
    y <- ifelse(x <= 0, 
                x * exp(delta_left/2 * x^2),
                x * exp(delta_right/2 * x^2))
}
lamw_untransform_hh <- function(x, delta_left, delta_right) {
    delta_left <- max(delta_left, 1e-16)
    delta_right <- max(delta_right, 1e-16)
    x <- ifelse(y <= 0, 
                sign(y) * sqrt(W(delta_left * y^2) / delta_left), 
                sign(y) * sqrt(W(delta_right * y^2) / delta_right))
}

x <- seq(-4, 4, 0.01)
y <- lamw_transform_h(x, delta)
plot(x,y,ylim=c(-40,40))
lines(x,x)
mtext('h transform, delta=1/3', side = 3, line = -1, adj = 0.025, padj = 1, cex = 1, col = "grey40")

y_inv <- lamw_untransform_h(y, delta)
all(abs(y_inv - x) < 1e-8)

x <- seq(-4, 4, 0.01)
y <- lamw_transform_hh(x, delta_left, delta_right)
plot(x,y,ylim=c(-40,40), yaxt='n')
lines(x,x)
mtext('hh transform, delta=(1/3, 0)', side = 3, line = -1, adj = 0.025, padj = 1, cex = 1, col = "grey40")

y_inv <- lamw_untransform_hh(y, delta_left, delta_right)
all(abs(y_inv - x) < 1e-8)

# Generated samples  ----------------------------------------------------------

quantile(as_draws_df(mod_out$draws())$new_y, c(0.05,0.25,0.5,0.75,0.95))
quantile(y, c(0.05,0.25,0.5,0.75,0.95))