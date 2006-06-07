data.restore("myeloma.dump")

myeloma <- myeloma[myeloma$entry < myeloma$futime,]

fit0 <- coxph(Surv(entry, futime, status) ~1, myeloma)
fit1 <- coxph(Surv(entry, futime, status) ~ dx.dt, myeloma)
fitx <- coxph(Surv(futime, status) ~ pspline(dx.dt, trace=T), myeloma)

# take a break while this one runs -- lunch maybe
fit2 <- coxph(Surv(entry, futime, status) ~ pspline(dx.dt, trace=T), myeloma)
