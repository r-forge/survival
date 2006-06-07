#
# Check out the dfbeta matrix from aareg
#   Note that it is kept internally in time order, not data set order
# Those who want residuals should use the resid function!

#
# First, the simple test case where I know the anwers
#
afit <- aareg(Surv(time, status) ~ x, test1, dfbeta=T)
temp <- c(rep(0,6),             #intercepts at time 1
          c(2,-1,-1,0,0,0)/9,   #alpha at time 1
          c(0,0,0,2, -1, -1)/9, #intercepts at time 2
          c(0,0,0,-2,1,1)/9)    #alpha at time 2
aeq(afit$dfbeta, temp)

#
#Now a multivariate data set
#
afit <- aareg(Surv(time, status) ~ age + sex + ph.ecog, lung, dfbeta=T)

ord <- order(lung$time, -lung$status)
cfit <- coxph(Surv(time, status) ~ age + sex + ph.ecog, lung[ord,],
	        method='breslow', iter=0, x=T)
cdt  <- coxph.detail(cfit, riskmat=T)

# an arbitrary list of times
acoef <- rowsum(afit$coef, afit$times) #per death time coefs
indx <- match(cdt$time, afit$times)
for (i in c(2,5,27,54,101, 135)) {
    lwho <- (cdt$riskmat[,i]==1)
    lmx <- cfit$x[lwho,]
    lmy <- 1*( cfit$y[lwho,2]==1 & cfit$y[lwho,1] == cdt$time[i])
    fit <- lm(lmy~ lmx)
    cat("i=", i, "coef=", aeq(fit$coef, acoef[i,]))

    rr <- diag(resid(fit))
    zz <- cbind(1,lmx)
    zzinv <- solve(t(zz) %*% zz)
    cat(" twt=", aeq(1/(diag(zzinv)), afit$tweight[indx[i],]))

    df <- t(zzinv %*% t(zz)  %*% rr)
    cat(" dfbeta=", aeq(df, afit$dfbeta[lwho,,i]), "\n")
    }
	  
rm(afit, cfit, cdt, lwho, lmx, lmy, fit, rr, zz, df)


# Repeat it with case weights
ww <- rep(1:5, length=nrow(lung))/ 3.0
afit <- aareg(Surv(time, status) ~ age + sex + ph.ecog, lung, dfbeta=T,
	      weights=ww)
cfit <- coxph(Surv(time, status) ~ age + sex + ph.ecog, lung[ord,],
	        method='breslow', iter=0, x=T, weight=ww[ord])
cdt  <- coxph.detail(cfit, riskmat=T)

acoef <- rowsum(afit$coef, afit$times) #per death time coefs
for (i in c(2,5,27,54,101, 135)) {
    who <- (cdt$riskmat[,i]==1)
    x <- cfit$x[who,]
    y <- 1*( cfit$y[who,2]==1 & cfit$y[who,1] == cdt$time[i])
    w <- cfit$weight[who]
    fit <- lm(y~x, weights=w)
    cat("i=", i, "coef=", aeq(fit$coef, acoef[i,]))

    rr <- diag(resid(fit))
    zz <- cbind(1,x)
    zzinv <- solve(t(zz)%*% (w*zz))
    cat(" twt=", aeq(1/(diag(zzinv)), afit$tweight[indx[i],]))
 
    df <- t(zzinv %*% t(zz) %*% (w*rr))
    cat(" dfbeta=", aeq(df, afit$dfbeta[who,,i]), "\n")
    }
	  
rm(afit, cfit, cdt, who, x, y, fit, rr, zz, df)
rm(ord, acoef)
    
