.BG
.FN survreg
.TL
Regression for a parametric survival model
.CS
survreg(formula, data=sys.parent(), subset, na.action,
link=c("log", "identity"),
dist=c("extreme", "logistic", "gaussian", "exponential"),
scale, eps=0.0001, init, iter.max=10, model=F, x=F, y=F, ...)
.AG formula
a formula expression as for other regression models.
See the documentation for `lm' and `formula' for details.
.AG data
optional data frame in which to interpret the variables occuring in the
formula.
.AG subset
subset of the observations to be used in the fit.
.AG na.action
function to be used to handle any NAs in the data.
.AG link
transformation to be used on the y variable.
.AG dist
assumed distribution for the transformed y variable.
.AG scale
optional; gives a fixed value for the scale parameter.
.AG eps
convergence criteria for the computation.  Iteration continues until the
relative change in log likelihood is less than eps.
.AG init
optional vector of initial values for the paramters.
.AG iter.max
maximum number of iterations to be performed.
.AG model
if TRUE, the model frame is returned.
.AG x
if TRUE, then the X matrix is returned.
.AG y
if TRUE, then the y vector (or survival times) is returned.
.AG ...
all the optional arguments to lm, including `singular.ok'.
.RT
an object of class `survreg' is returned, which inherits from class `glm'.
.SH Computation
  This routine is not as robust against nearly singular X matrices as lm();
the problem occurs when we explicitly invert the covariance matrix with
solve().  This can sometimes be solved by subtracting the mean from all
continuous covariates.
.EX
survreg(Surv(futime, fustat) ~ ecog.ps + rx, fleming, dist='extreme',
		link='log', scale=1)   #Fit an exponential
.KW survival
.WR
