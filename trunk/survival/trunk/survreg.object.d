.BG D
.FN survreg.object
.TL
Parametric Survival Model Object
.PP
This class of objects is returned by the `survreg' function
to represent a fitted parametric survival model.
Objects of this class have methods for the functions `print',
`summary', `predict', and 'residuals'.
.SH COMPONENTS
The following components must be included in a legitimate `survreg' object.
.AG coefficients
the coefficients of the `linear.predictors', which multiply  the
columns of the model
matrix.
It does not include the estimate of error (sigma).
The names of the coefficients are the names of the
single-degree-of-freedom effects (the columns of the
model matrix).
If the model is over-determined there will
be missing values in the coefficients corresponding to inestimable
coefficients.
.AG icoef
coefficients of the baseline model, which will contain the intercept
and log(scale), or mulitple scale factors for a stratified model.
.AG var
the variance-covariance matrix for the parameters, including the log(scale)
parameter(s).
.AG loglik
a vector of length 2, containing the log-likelihood for the baseline and
full models. 
.AG iter 
the number of iterations required
.AG linear.predictors
the linear predictor for each subject.
.AG df
the degrees of freedom for the final model.  For a penalized model
this will be a vector with one element per term.
.AG scale
the scale factor(s), with length equal to the number of strata.
.AG idf
degrees of freedom for the initial model.
.AG means
a vector of the column means of the coefficient matrix.
.AG dist
the distribution used in the fit.
.PP
The object will also have the following components found in 
other model results (some are optional):
`linear predictors', `weights', 'x', 'y', 'model', 
`call', `terms' and `formula'.
See `lm.object'.
.SA
`survreg'
.KW regression
.KW survival
.WR
