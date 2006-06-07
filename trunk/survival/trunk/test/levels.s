#
# Check out how predictions work when there are factors
#  in the model
#
fit1 <- coxph(Surv(time, status) ~ age + sex + factor(inst), lung)

tdata1 <- lung[1:30,]  #Notice that this is missing several levels
tdata2 <- tdata1; tdata2$inst <- factor(tdata2$inst)

# Note that pne subject is messing institution
imiss <- !is.na(lung$inst)
pred1 <- predict(fit1, type='lp', se=T)
pred2 <- predict(fit1, newdata=lung, type='lp', se=T)
pred3 <- predict(fit1, newdata=tdata1, type='lp', se=T)

aeq(pred1$fit[imiss], pred2$fit)
aeq(pred1$se[imiss], pred2$se)
aeq(pred1$fit[1:30], pred3$fit)
aeq(pred1$se[1:30],  pred3$se)


pred1 <- predict(fit1, type='terms', se=T)
pred2 <- predict(fit1, newdata=lung, type='terms', se=T)
pred3 <- predict(fit1, newdata=tdata1, type='terms', se=T)

aeq(pred1$fit[imiss,], pred2$fit)
aeq(pred1$fit[1:30,], pred3$fit)
