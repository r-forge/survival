.BG D
.FN coxph.object
.TL
Proportional Hazards Regression Object
.PP
This class of objects is returned by the `coxph' class of functions
to represent a fitted proportional hazards model.

Objects of this class have methods for the functions `print',
`summary', `residuals', `predict' and `survfit'.
.SH COMPONENTS
The following components must be included in a legitimate `coxph' object.
.AG coefficients
the coefficients of the linear predictor, which multiply the columns of the
model matrix.  If the model is over-determined there will be missing
values in the vector corresponding to the redundant columns in the model
matrix.
.AG var
the variance matrix of the coefficients.  Rows and columns corresponding to
any missing coefficients are set to zero.
.AG naive.var
this component will be present only if the `robust' option was true.  If so,
the `var' component will contain the robust estimate of variance, and this
component will contain the ordinary estimate.
.AG loglik
a vector of length 2 containing the log-likelihood with the initial values and
with the final values of the coefficients.
.AG score
value of the efficient score test, at the initial value of the coefficients.
.AG iter
number of iterations used.
.AG icc
the intra-class correlation coefficient.  This is only present if a
`cluster' term appeared in the model and the number of subjects per cluster
was >1.  It is a vector containing the number of clusters, the intra-class
correlation of the martingale residuals, and the icc of the ranked residuals.
Because the residuals are skewed, the latter may differ from the ordinary icc.
.AG linear.predictors
the vector of linear predictors, one per subject.
.AG residuals
the martingale residuals.
.AG means
vector of column means of the X matrix.  Subsequent survival curves are
adjusted to this value.
.AG n
the number of observations used in the fit.
.AG weights
the vector of case weights, if one was used.
.AG method
the computation method used.
.AG na.action
the na.action attribute, if any, that was returned by the na.action routine.
.PP
The object will also contain the following, for documentation see the `lm'
object: `terms', `assign', `formula', `call', and, optionally, `x', `y',
and/or `frame'.
.SA
`survfit', `coxph.detail', `cox.zph', `survreg', `residuals.coxph'.
.KW survival
.WR
