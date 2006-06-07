#
# Test aareg, for some simple data where the answers can be computed
#  in closed form
#  
test1 <- data.frame(time=  c(4, 3,1,1,2,2,3),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0),
		    wt=    c(1, 1:6))

tfit  <- aareg(Surv(time, status) ~ x, test1)
aeq(tfit$times, c(1,2,2))
aeq(tfit$nrisk, c(6,4,4))
aeq(tfit$coefficient, matrix(c(0,0,1/3, 1/3, 1, -1/3), ncol=2))
aeq(tfit$tweight, matrix(c(3,3,3, 3/2, 3/4, 3/4), ncol=2))
aeq(tfit$test.statistic, c(1,1))
aeq(tfit$test.var,       c(1, -1/4, -1/4, 1/4 + 9/16 + 1/16))

tfit <- aareg(Surv(time, status) ~ x, test1, test='nrisk')
aeq(tfit$tweight, matrix(c(3,3,3, 3/2, 3/4, 3/4), ncol=2)) #should be as before
aeq(tfit$test.statistic, c(4/3, 6/3+ 4 - 4/3))
aeq(tfit$test.var,       c(16/9, -16/9, -16/9, 36/9 + 16 + 16/9))

# In the 1-variable case, this is the same as the default Aalen weight
tfit <- aareg(Surv(time, status) ~ x, test1, test='variance')
aeq(tfit$test.statistic, c(1,1))
aeq(tfit$test.var,       c(1, -1/4, -1/4, 1/4 + 9/16 + 1/16))

#
# Repeat the above, with case weights
#
tfit <- aareg(Surv(time, status) ~x, test1, weights=wt)  
aeq(tfit$times, c(1,2,2))
aeq(tfit$nrisk, c(21,16,16))
aeq(tfit$coefficient, matrix(c(0,0,5/12, 2/9, 1, -5/12), ncol=2))
aeq(tfit$tweight, matrix(c(12,12,12, 36/7, 3,3), ncol=2))
aeq(tfit$test.statistic, c(5, 72/63 + 3 - 15/12))
aeq(tfit$test.var,       c(25, -25/4, -25/4, (72/63)^2 + 9 + (5/4)^2))

tfit <- aareg(Surv(time, status) ~x, test1, weights=wt, test='nrisk')
aeq(tfit$test.statistic, c(20/3,  42/9 + 16 - 16*5/12))
aeq(tfit$test.var,       c(400/9, -400/9, -400/9, 
			    (42/9)^2 + 16^2 + (16*5/12)^2))

#
# Make a test data set with no NAs, in sorted order, no ties,
#   15 observations
tdata <- lung[15:29, c('time', 'status', 'age', 'sex', 'ph.ecog')]
tdata$status <- tdata$status -1
tdata <- tdata[order(tdata$time, tdata$status),]
row.names(tdata) <- 1:15
tdata$status[8] <- 0      #for some variety

afit <- aareg(Surv(time, status) ~ age + sex + ph.ecog, tdata, nmin=6)
#
# Now, do it "by hand"
cfit <- coxph(Surv(time, status) ~ age + sex + ph.ecog, tdata, iter=0,
               method='breslow')
dt1   <- coxph.detail(cfit)
sch1  <- resid(cfit, type='schoen')

# First estimate of Aalen: from the Cox computations, first 9
#  The first and last cols of the ninth are somewhat unstable (approx =0)
mine <- rbind(solve(dt1$imat[,,1], sch1[1,]),
              solve(dt1$imat[,,2], sch1[2,]),
              solve(dt1$imat[,,3], sch1[3,]),
              solve(dt1$imat[,,4], sch1[4,]),
              solve(dt1$imat[,,5], sch1[5,]),
              solve(dt1$imat[,,6], sch1[6,]),
              solve(dt1$imat[,,7], sch1[7,]),
              solve(dt1$imat[,,8], sch1[8,]),
              solve(dt1$imat[,,9], sch1[9,])) 
mine <- diag(1/dt1$nrisk[1:9]) %*% mine

aeq(mine, afit$coef[1:9, -1])

rm(tfit, afit, mine, dt1, cfit, sch1)

