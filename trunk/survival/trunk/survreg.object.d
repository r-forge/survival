.BG D
.FN survreg.object
.TL
Parametric Survival Model Object
.PP
This class of objects is returned by the `survreg' function
to represent a fitted parametric survival model.
Class `survreg' inherits from class `glm', since it is fit by iterative
reweighted least squares; the object returned has all the components of a
weighted least squares object.

Objects of this class have methods for the functions `print',
`summary', `predict', and 'residuals'.
.SH COMPONENTS
The following components must be included in a legitimate `survreg' object.
The residuals, fitted values, coefficients and effects should be extracted
by the generic functions of the same name, rather than
by the `"$"' operator. 
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
.AG parms
the parameters of the model that are not coefficients of the X matrix.
The first of these will always be `log(sigma)'.
.AG fixed
a vector of the same length as `parms', where 1 indicates a parameter that
was fixed at its starting value and was not part of the iteration.
.AG deviance
minus twice the difference between the maximized log-likelihood under the
fitted model and a saturated model.
Similar to the residual sum of squares.
.AG loglik
the log-likelihood for the final model.
.AG null.deviance
the deviance corresponding to the model with only an itercept term, and
with `parms' fixed at their final values.
.AG dresiduals
the deviance residuals.
.AG var
the final variance matrix, including both coefficients and free parameters.
.AG family
a 2 element character vector giving the name of the family and
the link; mainly for printing purposes.
.PP
The object will also have the components of an `glm' object:
`linear predictors', `fitted.values', `residuals',
`effects', `R', `rank', `assign', `contrasts', `weights', `iter',
`residuals', `fitted.values', `call', `terms' and `formula'.
See `glm.object'.
.SA
`survreg', `glm.object', `lm.object'.
.KW regression
.KW survival
.WR
