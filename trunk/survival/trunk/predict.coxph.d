.BG
.FN predict.coxph
.TL
Make predictions from a fitted coxph object.
.DN
Creates predicted values from a coxph fit.
.CS
predict.coxph(object, newdata, type=c("lp", "risk", "expected", "terms"), 
se.fit=F, terms=labels.lm(object), collapse)
.RA
.AG object
the results of a coxph fit.
.OA
.AG newdata
data for prediction.  If absent predictions are for the
subjects used in the original fit.
.AG type
the type of predicted value.
Choices are the linear predictor (lp), the risk score exp(lp) (risk),
the expected number of events given the covariates and follow-up time
(expected), and the terms of the linear predictor (terms).
.AG se.fit
if TRUE, pointwise standard errors are produced for the predictions.
.AG terms
if type="terms", this argument can be used to specify which terms should be
included; the default is all.
.AG collapse
optional vector of subject identifiers. 
If specified, the output will contain one entry per subject rather than one
entry per observation.
.RT
a vector or matrix of predictions, or a list containing the predictions
(element "fit") and their standard errors (element "se.fit") if the se.fit
option is TRUE.
.DT
This function is a method for the generic function 'predict' for class
'coxph'.  It can be invoked by calling 'predict' for an object of the
appropriate class, or by calling predict.coxph directly.
.SA
coxph, residuals.coxph
.EX
> options(na.action='na.omit')
> fit <- coxph(Surv(time, status) ~ age + ph.ecog + strata(inst), lung)
> mresid <- lung$status - predict(fit, type='expected') #Martingale resid
.KW survival
.WR
