.BG
.FN residuals.coxph
.TL
Calculate residuals for a coxph fit.
.DN
Calculates martingale, deviance, score or schoenfeld residuals for
Cox or Anderson-Gill proportional hazards models.
.CS
residuals.coxph(object, 
       type=c("martingale", "deviance", "score", "schoenfeld"), 
       miss.expand=T, collapse)
.RA
.AG object
a coxph object, output from a coxph fit.
.OA
.AG type
character string indicating the type of residual desired.  Default is ???.
Only enough of the string to determine a unique match is required.
.AG miss.expand
Should the function return one residual for each row (miss.expand == T)
or for each row without NA's (miss.expend == F)?  Default is TRUE.
Note: argument is only relevant when data set contains missing values and 
when na.action=='naomit'.
.AG collapse
Vector indicating which rows to collapse(sum) over.  In time-dependent
models more than one row can pertain to the one individual.  If collapse
is given a vector indicating which rows correspond the the same individuals
one residual for each individual will be returned.
.RT
The object returned will be a vector for martingale and deviance 
residuals and matrices for score and schoenfeld residuals.  There will
be one row of residuals for each row in the input unless miss.expand=F
or collapse is specified.  One column of score and schoenfeld 
residuals will be returned for each column in the model.matrix.
.SE
~describe any side effects if they exist
.DT
~explain details here.
.SH REFERENCES
Terry Therneau, author of local function.

P. Andersen and R. Gill. "Cox's regression model for
counting processes, a large sample study", Annals of Statistics, 
10:1100-1120, 1982.  

T.Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.
.SA
coxph.
.EX

#
# An example of how to manage the residuals
#
attach(jasa1)
> fit <- coxph(Surv(start, stop, event) ~ (age + surgery)* transplant)
#
# The true residual for a subject is the sum of the rows that were used to
#   represent him/her to coxph.
> mresid <- resid(fit, collapse=jasa1$id)

.KW ~keyword
.WR
