.BG
.FN strata
.TL
Identify strata variables.
.DN
This is a special function used in the context of the Cox model.  It
identifies stratification variables when they appear on the right hand
side of a formula.
.CS
strata(..., na.group=F)
.RA
.AG ...
Any number of variables.  All must be the same length.
.OA
.AG na.group
if set to T, then missing values are treated as a distinct level of each
variable.
.RT
a new factor, whose levels are all possible combinations of the factors
supplied as arguments.
.DT
The result is identical to the interaction() function, but for the
labeling of the factors (strata is more verbose).
.SA
coxph
.EX
coxph(Surv(futime, fustat) ~ age + strata(rx))
.KW survival
.WR
