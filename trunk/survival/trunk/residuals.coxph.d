.BG
.FN residuals.coxph
.TL
Calculate residuals for a coxph fit.
.DN
Calculates martingale, deviance, score or Schoenfeld residuals for a
Cox proportional hazards model.
.CS
resid(object,
       type=c("martingale", "deviance", "score", "schoenfeld"), 
       collapse)
.RA
.AG object
a coxph object, output from a coxph fit.
.OA
.AG type
character string indicating the type of residual desired;
the default is martingale.
Only enough of the string to determine a unique match is required.
.AG collapse
Vector indicating which rows to collapse(sum) over.  In time-dependent
models more than one row data can pertain to a single individual.
If there were 4 individuals represented by 3, 1, 2 and 4 rows of data
respectively, then `collapse=c(1,1,1, 2, 3,3, 4,4,4,4)' could be used to
obtain per subject rather than per observation residuals.
.RT
The object returned will be a vector for martingale and deviance 
residuals and matrices for score and schoenfeld residuals.  There will
be one row of residuals for each row in the input data (without `collapse').
One column of score and Schoenfeld
residuals will be returned for each column in the model.matrix.
.SE
For deviance residuals, the status variable may need to be reconstructed.
For score and Shoenfeld residuals, the X matrix will need to be reconstructed.
.SH REFERENCES
T.Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.
.SA
coxph
.EX
> attach(jasa1)
> fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant)
> mresid <- resid(fit, collapse=jasa1$id)
.KW survival
.WR
