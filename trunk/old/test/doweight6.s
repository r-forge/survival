#
# Effect of weights on the robust variance
#
test1 <- data.frame(test1,
		    wt=    c(3,0,1,1,1,1,1))
testx <- data.frame(time=  c(4,4,4,1,1,2,2,3),
                    status=c(1,1,1,1,0,1,1,0),
                    x=     c(0,0,0,1,1,1,0,0),
		    wt=    c(1,1,1,1,1,1,1,1))
 
fit1 <- coxph(Surv(time, status) ~x, test1, method='breslow', weights=wt,
	      robust=T)
fit2 <- coxph(Surv(time, status) ~x, testx, method='breslow', robust=T)

db1 <- resid(fit1, 'dfbeta', weighted=F)
db1 <- db1[-2]         #toss the missing
db2 <- resid(fit2, 'dfbeta')
aeq(db1, db2[3:8])

W <- c(3,1,1,1,1,1)   #Weights, after removal of the missing value
aeq(fit2$var, sum(db1*db1*W))
aeq(fit1$var, sum(db1*db1*W*W))

