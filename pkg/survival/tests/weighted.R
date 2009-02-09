# Tests of the weighted Cox model, AG form of the data
#   Same solution as book5 and book6
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

testw3 <- data.frame(id  =  c( 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 8, 8, 9),
		     begin= c( 0, 5, 0, 0,10,15, 0, 0,14, 0, 0, 0,23, 0),
		     time=  c( 5,10,10,10,15,20,20,14,20,20,30,23,40,50),
		    status= c( 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0),
		    x=      c( 2, 2, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0),
		    wt =    c( 1, 1, 2, 3, 3, 3, 4, 3, 3, 2, 1, 2, 2, 1))

fit0 <- coxph(Surv(begin,time, status) ~x, testw3, weights=wt,
		    method='breslow', iter=0)
fit  <- coxph(Surv(begin,time, status) ~x, testw3, weights=wt, method='breslow')
fit0
summary(fit)
resid(fit0, type='mart', collapse=testw3$id)
resid(fit0, type='score', collapse=testw3$id)
resid(fit0, type='scho')

resid(fit, type='mart', collapse=testw3$id)
resid(fit, type='score', collapse=testw3$id)
resid(fit, type='scho')
fit0 <- coxph(Surv(begin, time, status) ~x,testw3, weights=wt, iter=0)
resid(fit0, 'mart', collapse=testw3$id)
resid(coxph(Surv(begin, time, status) ~1, testw3, weights=wt)
		      , collapse=testw3$id)  #Null model

fit  <- coxph(Surv(begin,time, status) ~x, testw3, weights=wt, method='efron')
fit
resid(fit, type='mart', collapse=testw3$id)
resid(fit, type='score', collapse=testw3$id)
resid(fit, type='scho')

#
# Effect of case weights on expected survival curves post Cox model
#
fit0  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='breslow',
	       iter=0)
fit0b <- coxph(Surv(time, status) ~x, testw2, method='breslow', iter=0)

surv1 <- survfit(fit0, newdata=list(x=0))
surv2 <- survfit(fit0b, newdata=list(x=0))
aeq(surv1$surv, surv2$surv)


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

