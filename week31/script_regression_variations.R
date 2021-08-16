library(cmdstanr)
library(posterior)
library(LambertW)

# setup -----------------------------------------------------------------------

N <- 1000

lambda <- 3
alpha <- 1
beta <- 3
x <- rpois(N, lambda=lambda)

mu <- beta*x + a
sigma <- 3/2

error1 <- function(sigma) {
    rnorm(N, 0, sigma)
}

error2 <- function(sigma) {
    nu <- 2*sigma^2/(sigma^2 - 1)
    rt(N, nu)
}

error3 <- function(sigma, delta=1/3) {
    u <- rnorm(N)
    u*exp(delta/2*u^2)*sigma
}

error4 <- function(sigma) {
    e1 <- rnorm(N, 0, sigma)
    nu <- 2*sigma^2/(sigma^2 - 1)
    e2 <- rt(N, nu)
    ifelse(e1 >= 0, e1, -abs(e2))
}

e1 <- error1(sigma)
e2 <- error2(sigma)
e3 <- error3(sigma)
e4 <- error4(sigma)
e <- c(e1, e2, e3, e4)
dim(e) <- c(N, 4)

# fit stan program ------------------------------------------------------------

fp <- file.path(paste(getwd(), "/week31/regression_lambertw_normal_hh.stan", sep=""))
mod <- cmdstan_model(fp, force_recompile = F)

for (i in 1:4) {
    y <- alpha + beta*x + e[,i]
    
    mod_out <- mod$sample(data=list(N=N, y=y, x=x), parallel_chains=4)
    print(mod_out$summary()[1:6,])
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