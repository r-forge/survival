.BG
.FN survexp
.TL
Compute Expected Survival
.DN
Returns either the expected survival of a cohort of subjects, or the
individual expected survival for each subject.
.CS
survexp(formula, data, weights, subset, na.action,
 times, cohort=T, conditional=F,
 ratetable=survexp.us, scale=1, se.fit, model=F, x=F, y=F)
.RA
.AG formula
a formula object.  The response variable will be a vector of follow-up times,
and is optional.  The predictors will consist of optional grouping variables
separated by + operators (exactly as in `survfit'), along with a
`ratetable()' term.  This latter matches each subject to his/her expected
cohort.
.OA
.AG data, weights, subset, na.action
as in other modeling routines.
Weights are currently ignored.
.AG times
an optional vector of times at which the resulting survival curve should
be evaluated.  If absent, the result will be reported for each unique value
of the vector of follow-up times.
.AG cohort
If false, each subject is treated as a subgroup of size 1.
.AG conditional
If `y' is missing in the formula, this argument is ignored.  Otherwise it
is an indicator of whether y includes death times, which leads to conditional
expected survival, or y includes only the potential censoring times.
.AG ratetable
a table of event rates, such as survexp.uswhite, or a fitted Cox model.
.AG scale
a scaling for the results.  As most rate tables are in units/day, a
value of 365.24 would cause the output to be reported in years.
.AG npoints
calculate intermediate results at npoints values, evenly spaced on the range
of `y'.  The usual (exact) calculation is done at each unique 'y' value;
for very large data sets this may incur too much storage for the scratch
array.
For a prediction from a Cox model this arument is ignored.
.AG se.fit
compute the standard error of the predicted survival.
The default is to compute this whenever the routine can, which at this time
is only for the Ederer method and a Cox model as the rate table.
.AG model, x, y
flags to control what is returned.  If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result, with the same names as the
flag arguments.
.RT
if cohort=T an object of class `survexp', otherwise a vector of per-subject
expected survival values.  The former contains the number of subjects at
risk and the expected survival for the cohort at each requested time.
.DT
Individual expected survival is ususally used in models or testing, to
`correct' for the age and sex composition of a group of subjects.  For
instance, assume that birth date, entry date onto the study,sex and
actual survival time are all known for a group of subjects.
The uswhite population tables contain expected death rates
based on calendar year, sex and age.  Then
.Cs
haz <- -log(survexp(death.time ~ ratetable(sex=sex, year=entry.dt, age=(birth.dt-entry.dt)), cohort=F))
.Ce
gives for each subject the total hazard experienced up to their observed
death time or censoring time.
This probability can be used as a rescaled time value in models:
.Cs
glm(status ~ 1 + offset(log(haz)), family=poisson)
glm(status ~ x + offset(log(haz)), family=poisson)
.Ce
In the first model, a test for intercept=0 is the one sample log-rank
test of whether the observed group of subjects has equivalent survival to
the baseline population.  The second model tests for an effect of variable
`x' after adjustment for age and sex.
.PP
Cohort survival is used to produce an overall survival curve.  This is then
added to the Kaplan-Meier plot of the study group for visual comparison
between these subjects and the population at large.  There are three common
methods of computing cohort survival.
In the "exact method" of Ederer the cohort is not censored; this corresponds
to having no response variable in the formula.  Hakulinen recommends censoring
the cohort at the anticipated censoring time of each patient, and Verhuel
recommends censoring the cohort at the actual observation time of each
patient.
The last of these is the conditional method.
These are obtained by using the respective time values as the
follow-up time or response in the formula.
.SH REFERENCES
G. Berry.  The analysis of mortality by the subject-years method.
Biometrics 1983, 39:173-84.
.br
F Ederer, L Axtell, and S Cutler.  The relative survival rate: a statistical
methodology. Natl Cnacer Inst Monogr 1961, 6:101-21.
.br
T. Hakulinen.  Cancer survival corrected for heterogeneity in patient
withdrawal.  Biometrics 1892, 38:933.
.br
H. Verheul, E. Dekker, P. Bossuyt, A. Moulijn, and A. Dunning.  Backround
mortality in clinical survival studies.  Lancet 1993, 341:872-5.
.SA
survfit, survexp.us, survexp.fit, personyr, date
.EX
efit <- survexp( ~ ratetable(sex=sex, year=entry.dt, age=entry.dt-birth.dt))
plot(survfit(Surv(futime, status) ~1))
lines(efit)
.KW survival
.WR
