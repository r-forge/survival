# This data set is the result of a simulation, done by
#   Jose Cortinas at EORTC.
# It has a random center effect, and a random treatment effect
load('eortc.rda')

fit0 <- coxph(Surv(y, uncens)~ trt, eortc)  #Simple fit
fit1 <- coxme(Surv(y, uncens) ~ trt + (1|center), eortc)

# Random treatment effect, + center effect
#   fit2b by default assumes some sparseness, a and c do not
fit2a <- coxme(Surv(y, uncens) ~trt + (1| center/trt), eortc, sparse=c(1,0))
fit2b <- coxme(Surv(y, uncens) ~trt + (1| center/trt), eortc)
fit2c <- coxme(Surv(y, uncens) ~trt + (1| center/trt), eortc,
               varlist=coxvarFull(collapse=TRUE))
all.equal(fit2a$loglik, fit2b$loglik, tolerance=1e-5)
all.equal(fit2c$loglik, fit2a$loglik)
 aeq(fixef(fit2a), fixef(fit2c))

# Treatment as a random slope
fit3 <- coxme(Surv(y, uncens) ~ trt + (1| center) + (trt|center), eortc)

# Same model as fit2c, using matrices
group <- paste(eortc$center, eortc$trt, sep='/')
ugroup <- paste(rep(1:37, each=2), rep(0:1, 37), sep='/') #unique groups
mat1 <- bdsmatrix(rep(c(1,1,1,1), 37), blocksize=rep(2,37),
                  dimnames=list(ugroup,ugroup))
mat2 <- as.matrix(bdsI(ugroup))

fit4 <- coxme(Surv(y, uncens) ~trt + (1|group), eortc,
              varlist=coxvarMlist(list(mat1, mat2), rescale=F, pdcheck=F),
               vinit=c(.2, .2))
aeq(fit4$log, fit2c$log)
aeq(as.matrix(fit4$var), as.matrix(fit2c$var))
fit3 <- coxme(Surv(y, uncens) ~ trt + (1| center) + (trt|center), eortc)
fit3b <- coxme(Surv(y, uncens) ~ trt + (1| center/trt), eortc)


#
# Random treatment effect, correlated with the random intercept
#
mat3 <-  bdsmatrix(rep(c(0,1,1,1), 37), blocksize=rep(2,37),
                  dimnames=list(ugroup,ugroup))
fit3 <- coxme(Surv(y, uncens) ~x, eortc,
              random= ~1|group, varlist=list(mat1, mat2, mat3),
              rescale=F, pdcheck=F, vinit=c(.04, .12, .02))


# Random treatment, done differently
fit2b <- coxme(Surv(y, uncens) ~x, eortc,
              random= ~1|center/x, varlist=list(bdsI, bdsI))

mat2b <- bdsI(ugroup, blocksize=74)  #force it to not be sparse
fit2c <- coxme(Surv(y, uncens) ~x, eortc, 
              random= ~1|group, varlist=list(mat1, mat2b))

temp <- fit2c$coef$random
fit3c <- coxme(Surv(y, uncens) ~x, eortc,
              random= ~1|group, varlist=list(mat1, mat2, mat3),
              rescale=F, pdcheck=F, lower=c(0,0,-100),
              variance=c(temp[1]+temp[2], 2*temp[2], -temp[2]))


