.BG
.FN survdiff
.TL
Test Survival Curve Differences
.DN
Tests if there is a difference between two or more survival curves using
the G-rho family of tests, or for a single curve against a known alternative.
.CS
survdiff(formula, data, subset, na.action, rho=0)
.RA
.AG formula
a formula expression as for other survival models, of the form
`Surv(time, status) ~ predictors'.  For a one-sample test, the predictors
must consist of a single `offset(sp)' term, where `sp' is a vector giving the
survival probability of each subject.  For a k-sample test, each unique
combination of predictors defines a subgroup.
A `strata' term may be used to produce a stratified test.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the `strata' function with its
`na.group=T' argument.
.OA
.AG data
an optional data frame in which to interpret the variables occurring in the
formula.
.AG subset
expression indicating which subset of the rows of data should be used in
the fit.  This can be a logical vector (which is replicated to have
length equal to the number of observations), a numeric vector indicating
which observation numbers are to be included (or excluded if negative),
or a character vector of row names to be included.  All observations are
included by default.
.AG na.action
a missing-data filter function.  This is applied to the `model.frame' after any
subset argument has been used.  Default is `options()$na.action'.
.AG rho
a scalar parameter that controls the type of test.
.RT
a list with components:
.AG n
the number of subjects in each group.
.AG obs
the weighted observed number of events in each group.
If there are strata, this will be a matrix with one column per stratum.
.AG exp
the weighted expected number of events in each group.
If there are strata, this will be a matrix with one column per stratum.
.AG chisq
the chisquare statistic for a test of equality.
.AG var
the variance matrix of the test.
.AG strata
optionally, the number of subjects contained in each stratum.
.SH METHOD
This function implements the G-rho family of
Harrington and Fleming (1982), with weights on each death of `(S(t))^rho',
where `S' is the Kaplan-Meier estimate of survival.
When `rho = 0' this is the log-rank or Mantel-Haenszel test,
and when `rho = 1' it is equivalent to the Peto & Peto modification
of the Gehan-Wilcoxon test.
.PP
If the right hand side of the formula consists only of an offset term,
then a one sample test is done.
To cause missing values in the predictors to be treated as a separate
group, rather than being omitted, use the `factor' function with its
`exclude' argument.
.SH REFERENCE
Harrington, D. P. and Fleming, T. R. (1982).
A class of rank test procedures for censored survival data.
.ul
Biometrika
\fB69\fR, 553-566.
.EX
survdiff(Surv(futime, fustat) ~ rx)
survdiff(Surv(time, status) ~ pat.karno + strata(inst), data=cancer)

expect <- survexp(entry, birth, sex, futime)
survdiff(Surv(futime, fustat) ~ offset(expect$surv))  #One sample log-rank
.KW survival
.WR
