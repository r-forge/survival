.BG
.FN frailty
.TL
fit a penalized factor variable
.DN
fit a penalized factor variable in a survival model; random effects
or 'frailties' are one instance of this
.CS
frailty(x, distribution="gamma", ...)
.RA
.AG x
a variable describing discrete groups
.OA 
.AG distribution
the distribution of the random effect.
Currently, the legal values are gamma, gaussian, and t.
.AG ...
other arguments to the specific frailty function
.RT
a variable with attached attributes, so as to be recognized
as a penalized term by coxph or survreg
.DT
this routine calls frailty.gamma, frailty.gaussian, or frailty.t to do
the actual work.
.SA
coxph, survreg, frailty.gamma, frailty.gaussian, frailty.t
.EX
coxph(Surv(time, status) ~ age + sex + frailty(inst), lung
.KW survival
.WR
