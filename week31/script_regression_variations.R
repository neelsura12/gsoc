library(cmdstanr)
library(posterior)

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

# compare generated samples against originals

quantile(as_draws_df(mod_out$draws())$new_y, c(0.05,0.25,0.5,0.75,0.95))
quantile(y, c(0.05,0.25,0.5,0.75,0.95))