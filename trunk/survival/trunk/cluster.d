.BG
.FN cluster
.TL
Identify clusters.
.DN
This is a special function used in the context of the Cox model.  It
identifies correlated groups of observations, and is used on the right hand
side of a formula.
.CS
cluster(x)
.RA
.AG x
A character, factor, or numeric variable.
.RT
x
.DT
The function's only action is semantic, to mark a variable as the
cluster indicator.
.SA
coxph
.EX
coxph(Surv(futime, fustat) ~ age + cluster(group))
.KW survival
.WR
