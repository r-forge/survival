.BG
.FN survreg
.TL
Regression for a Parametric Survival Model
.DN 
Regression for a parametric survival model
.CS
survreg(formula, data=sys.parent(), subset, na.action,
dist="weibull", init, scale=0, control,
model=F, x=F, y=T, ...)
.AG formula
a formula expression as for other regression models.
See the documentation for `lm' and `formula' for details.
.OA
.AG data
optional data frame in which to interpret the variables occurring in the
formula.
.AG weights
case weights for the observations.
.AG subset
subset of the observations to be used in the fit.
.AG na.action
function to be used to handle any NAs in the data.
.AG dist
assumed distribution for y variable.  These include
weibull, exponential, gaussian, logistic, lognormal and loglogistic.
If the argument is a character string, then it is assumed to name an
element from `survreg.distributions' (enter `names(survreg.distributions)'
for a full list of choices).
Otherwise, it is assumed to be a user defined list conforming to this
standard.
.AG parm
a list of fixed parameters.  For the t-distribution for instance this is
the degrees of freedom; most of the distributions have no parameters.
.AG init
optional vector of initial values for the parameters.
.AG scale
optional fixed value for the scale.  If set to <=0 then the scale is
estimated.
.AG control
a list of control values, in the format producted by `survreg.control'.
.AG model
if `TRUE', the model frame is returned.
.AG x
if `TRUE', then the X matrix is returned.
.AG y
if `TRUE', then the y vector (or survival times) is returned.
.AG ...
other arguments which will be passed to `survreg.control'.
.RT
an object of class `survreg' is returned.
.SH NOTE
This routine underwent significant changes from survival4 to survival5.
.SH COMPUTATION
The routine uses a Newton-Raphson iteration with step halving, 
with provision for general penalized term.
Fisher scoring is used for intermediate steps where the information matrix
is not positive definite.
.SA
`survreg.object'
.EX
survreg(Surv(futime, fustat) ~ ecog.ps + rx, fleming, dist='lognormal')
.KW survival
.WR
