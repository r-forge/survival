#
# Check out the survreg density and probability functions
#

# Gaussian
x <- -10:10
p <- seq(.1, .95, length=25)
all.equal(dsurvreg(x, 1, 5, 'gaussian'), dnorm(x, 1, 5))
all.equal(psurvreg(x, 1, 5, 'gaussian'), pnorm(x, 1, 5))
all.equal(qsurvreg(p, 1, 5, 'gaussian'), qnorm(p, 1, 5))

# Lognormal
x <- 1:10
all.equal(dsurvreg(x, 1, 5, 'lognormal'), dlnorm(x, 1, 5))
all.equal(psurvreg(x, 1, 5, 'lognormal'), plnorm(x, 1, 5))
all.equal(qsurvreg(p, 1, 5, 'lognormal'), qlnorm(p, 1, 5))

# Weibull
lambda <- exp(-2)
rho    <- 1/3
temp <- (lambda*x)^rho
all.equal(psurvreg(x, 2, 3), 1- exp(-temp))
all.equal(dsurvreg(x, 2, 3), lambda*rho*(lambda*x)^(rho-1)*exp(-temp))
