.BG
.FN surv.diff
.TL
Test Survival Curve Differences
.DN
Tests if there is a difference between two or more survival curves using
the G-rho family of tests, or for a single curve against a known alternative.
.CS
surv.diff(formula, rho=0, subset)
.RA
.AG formula
a formula expression as for other survival models, of the form
`Surv(time, status) ~ predictors'.  For a one-sample test, the predictors
must consist of a single `offset(sp)' term, where sp is a vector giving the
survival probability of each subject.  For a k-sample test, each unique
combination of predictors defines a subgroup.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the `strata' function with its
`na.group=T' argument.
.OA
.AG rho
a parameter that controls the type of test.
.AG subset
subset of the observations to be used in the fit.
.RT
a list with components:
.AG n
the number of subjects in each group.
.AG obs
the weighted observed number of events in each group.
.AG exp
the weighted expected number of events in each group.
.AG chisq
the chisquare statistic for a test of equality.
.SH METHOD
This function implements the G-rho family of
Harrington and Fleming (1982), with weights on each death of (S(t))^rho,
where S is the Kaplan-Meier estimate of survival.
When `rho = 0' this is the log-rank or Mantel-Haenszel test,
and when `rho = 1' it is equivalent to the Peto & Peto modification
of the Gehan-Wilcoxon test.
.SH REFERENCE
Harrington, D. P. and Fleming, T. R. (1982).
A class of rank test procedures for censored survival data.
.ul
Biometrika
\fB69\fR, 553-566.
.SA
`surv.diff.print'.
.EX
surv.diff(Surv(futime, fustat) ~ rx)

expect <- surv.exp(entry, birth, sex, futime)
surv.diff(Surv(futime, fustat) ~ offset(expect$surv))  #One sample log-rank
.KW survival
.WR
