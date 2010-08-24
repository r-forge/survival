#
# Check out the code on a simulated data set.  On first cut it appeared to
#  do very badly, this turned out to be a bug in coxfit4d.c, when there
#  were random slopes it miscalculated the linear predictor.
#
library(coxme)
set.seed(1978)
mkdata <- function(n, beta=c(.4, .1), sitehaz=c(.5,1.5, 2,1)) {
    nsite <- length(sitehaz)
    site <- rep(1:nsite, each=n)
    trt1 <- rep(0:1, length=n*nsite)
    hazard <- sitehaz[site] + beta[1]*trt1 + beta[2]*trt1 * (site-mean(site))
    stime <- rexp(n*nsite, exp(hazard))
    q80 <- quantile(stime, .8)
    data.frame(site=site,
               trt1 = trt1,
               trt2 = 1-trt1,
               futime= pmin(stime, q80),
               status= ifelse(stime>q80, 0, 1),
               hazard=hazard
               )
    }

trdata <- mkdata(100)

set.seed(50)
nsim <- 500
fit <- coxme(Surv(futime, status) ~ trt2 + (1 + trt2 | site), trdata,
             refine.n=nsim)
debug <- fit$refine.debug

# Create the variance-covariance matrix that was used
set.seed(50)
bmat <- matrix(rnorm(8*nsim), nrow=8)
if (!is.null(debug)) all.equal(bmat, debug$b)

temp <- ranef(fit)[[1]]
sigma <- diag(c(rep(temp[1],4), rep(temp[4],4)))
sigma[cbind(1:4,5:8)] <- temp[3]* sqrt(temp[1] * temp[4])
sigma[cbind(5:8,1:4)] <- temp[3]* sqrt(temp[1] * temp[4])

if (!is.null(debug))
    all.equal(as.matrix(gchol(sigma), ones=F), as.matrix(debug$gkmat, ones=F))

b2 <- gchol(sigma) %*% bmat
coxll <- double(nsim)
for (i in 1:nsim) {
    lp <- trdata$trt2 * fixef(fit) + b2[trdata$site, i] +
           b2[trdata$site +4, i] * trdata$trt2
    tfit <- coxph(Surv(futime, status) ~ offset(lp), trdata)
    coxll[i] <- tfit$loglik
    }
if (!is.null(debug)) all.equal(coxll, debug$log)

# How does the direct average compare to the IPL?
constant <- .5*(log(2*pi) + sum(log(diag(gchol(sigma)))))
fit$log[2] + c(IPL=0, sim=log(mean(exp(coxll-fit$log[2]))) - constant)

# The better estimate uses control variables
bhat <- unlist(fit$frail)
temp <- t(b2-bhat) %*% fit$hmat[1:8,1:8]
laplace <-fit$log[3] -  .5*(rowSums(temp^2) - colSums(bmat^2))
#plot(coxll, laplace); abline(0,1)

#  To match the internal calculations of coxme, I need the uncorrected
#  IPL as a centering constant
fit2 <-  coxme(Surv(futime, status) ~ trt2 + (1 + trt2 | site), trdata)
errhat <- exp(coxll-fit2$log[2]) - exp(laplace-fit2$log[2])

all.equal(c(correction=mean(errhat), std=sqrt(var(errhat)/nsim)), fit$refine)
