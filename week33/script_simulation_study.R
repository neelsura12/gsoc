library(cmdstanr)
library(posterior)
library(LambertW)

# setup -----------------------------------------------------------------------

lamw_error <- function(N, sigma, delta_left, delta_right) {
    u <- rnorm(N)
    ifelse(u <= 0, u*exp(delta_left/2*u^2)*sigma, 
                   u*exp(delta_right/2*u^2)*sigma)
}

N <- 1000

mu_x <- 3
x <- rnorm(N, mu_x, 1)

alpha <- 1
beta <- 3
sigma <- 3/2

y1 <- alpha + beta*x + lamw_error(N, sigma, 0, 1/3)
y2 <- alpha + beta*x + lamw_error(N, sigma, 2/3, 1/3)

y <- c(y1, y2)
dim(y) <- c(N, 2)

# fit stan program ------------------------------------------------------------

fp <- file.path(paste(getwd(), "/week33/regression_lambertw_normal_hh.stan", sep=""))
mod <- cmdstan_model(fp, force_recompile = F)

mod_out = c()
for (i in 1:2) {
    tmp <- mod$sample(data=list(N=N, y=y[,i], x=x), parallel_chains=4)
    mod_out <- c(mod_out, tmp)
}

# make plots ------------------------------------------------------------------
# error histograms

par(mfrow = c(2, 2))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for (i in 1:4) {
    hist(e[,i], axes = T, freq = F, main=c(), yaxt="n")
    padj <- ifelse(i %in% c(3, 4), 1.5, 1)
    mtext(i, side = 3, line = -1, adj = 0.025, padj = padj, cex = 2, col = "grey40")
    box(col = "grey60")
}

# error qqplots

par(mfrow = c(2, 2))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
for (i in 1:4) {
    qqnorm(e[,i], pch = 1, frame = FALSE, main=c())
    qqline(e[,i], col = "steelblue", lwd = 2)
    padj <- ifelse(i %in% c(3, 4), 1.5, 1)
    mtext(i, side = 3, line = -1, adj = 0.025, padj = padj, cex = 2, col = "grey40")
    box(col = "grey60")
}

# LambertW  ------------------------------------------------------------------

par(mfrow = c(1, 2))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

hist(rnorm(N, 0, 1), main=c())
lamw <- rLambertW(N, "normal", theta=list(beta=c(0,1), delta=c(0,1/3), gamma=0, alpha=1))
hist(lamw, main=c())

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