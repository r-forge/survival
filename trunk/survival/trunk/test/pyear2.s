#
# Test out the case weights in pyears
#  We purposely force some of the data into the offtable slot
#
tdata <- data.frame(age=1:100 * 365, 
                    year=julian(sample(1:12, 100, replace=T),
                                rep(15,100),
                                sample(1930:2000, 100, replace=T)),
                    sex = rep(1:2, 50),
                    race= rep(1:2, length=100),
                    group = factor(rep(1:4,25)),
                    futime= sample(200:500, 100))

twt <- sample(1:6, 100, replace=T)
temp.age <- tcut(tdata$age, 20:55 * 365.25)
fit1 <- pyears(futime ~ group + temp.age, ratetable=survexp.us, data=tdata,
                 weight=twt)
repcount <-rep(1:100, twt) 
fit2 <- pyears(futime ~ group + temp.age, ratetable=survexp.us,
               data=tdata, subset=repcount)

aeq(fit1$pyears,   fit2$pyears)
# aeq(fit1$n,        fit2$n)   These aren't supposed to match
aeq(fit1$offtable, fit2$offtable)
aeq(fit1$expected, fit2$expected)


# Repeat, without expected (uses a different C routine)
fit1 <- pyears(futime ~ group + temp.age, data=tdata,  weight=twt)
fit2 <- pyears(futime ~ group + temp.age, data=tdata,   
               subset=rep(1:100, twt))

aeq(fit1$pyears,   fit2$pyears)
aeq(fit1$offtable, fit2$offtable)
