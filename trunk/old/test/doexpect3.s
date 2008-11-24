#
# A reported problem with a larger data set
#
tdata <- data.frame(age=rep(11:35,4) * 100*pi,
                    year = mdy.date(1,1,1950) + 1:100 * 200,
                    sex= rep(c('male', 'female'), 50),
                    id=sample(1:1000, 100))
tdata$birth.dt <- tdata$year - floor(tdata$age)

exp1 <- survexp( ~ id, data=tdata, ratetable=survexp.wnc,
                time= 1:20 * 182)$surv
exp2 <- exp1
idlist <- sort(tdata$id)
for (i in 1:100) {
    exp2[,i] <- survexp( ~ 1, data=tdata, ratetable=survexp.wnc,
                        subset=(id==idlist[i]), time=1:20 * 182)$surv
    }
aeq(exp1, exp2)


# But it does not work with this data set!
#
attach("../../../cuminc")
tdata <- mm2

tdata$age <- round(tdata$age*365.25)
tdata$sex <-  ifelse(tdata$sex=='M', 'male', 'female')
tdata$year <- tdata$roch.dt

fit1 <- survexp( ~clinic, tdata, ratetable=survexp.wnc,
                times= c(1:364, 365*1:10))

clist <- sort(unique(tdata$clinic))
fit2 <- survexp( ~1, tdata, ratetable=survexp.wnc,
                times=c(1:364, 365*1:10), subset=(clinic==clist[3]))
                 
aeq(fit1$surv[,3], fit2$surv)
haz1 <- diff(c(0, -log(fit1$surv[,3])))
haz2 <- diff(c(0, -log(fit2$surv)))
