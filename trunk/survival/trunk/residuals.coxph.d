.BG
.FN residuals.coxph
.TL
Calculate Residuals for a `coxph' Fit
.DN
Calculates martingale, deviance, score or Schoenfeld residuals for a
Cox proportional hazards model.
.CS
resid(object,
       type=c("martingale", "deviance", "score", "schoenfeld",
	      "dfbeta", "dfbetas", "scaledsch"),
       collapse=F, weighted=F)
.RA
.AG object
an object inheriting from class `coxph', representing a fitted Cox
regression model.
Typically this is the output from the `coxph' function.
.OA
.AG type
character string indicating the type of residual desired.
Possible values are `"martingale"', `"deviance"', `"score"', `"schoenfeld"',
"dfbeta"', `"dfbetas"', and `"scaledsch"'.
Only enough of the string to determine a unique match is required.
.AG collapse
vector indicating which rows to collapse (sum) over.
In time-dependent models more than one row data can pertain
to a single individual.
If there were 4 individuals represented by 3, 1, 2 and 4 rows of data
respectively, then `collapse=c(1,1,1, 2, 3,3, 4,4,4,4)' could be used to
obtain per subject rather than per observation residuals.
.AG weighted
if `TRUE' and the model was fit with case weights, then the weighted
residuals are returned.  
The default is to return unweighted residuals, except for dfbeta where
the default is to return weighted residuals.
.RT
For martingale and deviance residuals, the returned object is a vector
with one element for each subject (without `collapse').
For score residuals it is a matrix
with one row per subject and one column per variable.
The row order will match the input data for the original fit.
For Schoenfeld residuals, the returned object is a matrix with one row
for each event and one column per variable.  The rows are ordered by time
within strata, and an attribute `strata' is attached that contains the
number of observations in each strata.
The scaled Schoenfeld residuals are used in the `cox.zph' function.
.PP
The score residuals are each individual's contribution to the score vector.
Two transformations of
this are often more useful: `dfbeta' is the approximate change in the
coefficient vector if that observation were dropped,
and `dfbetas' is the approximate change in the coefficients, scaled by
the standard error for the coefficients.
.SH NOTE
For deviance residuals, the status variable may need to be reconstructed.
For score and Schoenfeld residuals, the X matrix will need to be reconstructed.
.SH REFERENCES
T. Therneau, P. Grambsch, and T. Fleming. "Martingale based residuals
for survival models", \fIBiometrika\fP, March 1990.
.SA
`coxph'
.EX
> attach(jasa1)
> fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant,
               data=heart)
> mresid <- resid(fit, collapse=heart$id)
.KW survival
.WR
