#
# Test out the survfit.ci function, which does competing risk
#   estimates
#
# For this we need the MGUS data set
#
data.restore('mgus1.sdump')

endpoint <- ifelse(is.na(mgus1$sdtime), 'death', 'progression')
fit1 <- survfit.ci(Surv(time, status) ~ strata(endpoint), mgus1)

# Now get the overall survival, and the hazard for progression
fit2 <- survfit(Surv(time, status) ~1, mgus1)
fit3 <- survfit(Surv(time, status*(endpoint=='progression')) ~1, mgus1,
                type='fleming')

# Classic CI formula
haz <- diff(c(0, -log(fit3$surv))) #Aalen hazard estimate
tsurv <- c(1, fit2$surv[-length(fit2$surv)])  #lagged survival
ci <- cumsum(haz *tsurv)
aeq(1-ci, fit1$surv[,2])

#
# Now, make sure that it works for subgroups
#
fit1 <- survfit.ci(Surv(time, status) ~ sex + strata(endpoint), mgus1)
fit2 <- survfit.ci(Surv(time, status) ~ 1 + strata(endpoint), mgus1,
                        subset=(sex==1))
fit3 <- survfit.ci(Surv(time, status) ~ 1 + strata(endpoint), mgus1,
                   subset=(sex==2))

aeq(fit2$surv, fit1$surv[1:fit1$strata[1],])
aeq(fit3$surv, fit1$surv[-(1:fit1$strata[1]),])

rm(fit1, fit2, fit3, tsurv, ci, haz, endpoint)
