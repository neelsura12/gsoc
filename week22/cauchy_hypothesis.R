set.seed(23)

dat <- rcauchy(100)
n <- length(dat)

x0 <- median(dat)
gamma <- IQR(dat)
tstat <- x0/(gamma/sqrt(n))

dat_normzld <- Gaussianize(dat)
tstat_normlzd <- mean(dat_normzld)/(sd(dat_normzld)/sqrt(n))
