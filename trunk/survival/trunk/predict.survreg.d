.BG
.FN predict.survreg
.TL
Predicted Values for a `survreg' Object
.DN
Predicted values for a `survreg' object
.CS
predict.survreg(object, newdata, 
type=c("response", "link", "linear", "response", "terms", "quantile", 
	"uquantile"), 
se.fit=F, terms=labels.lm(object), p=c(0.1, 0.9), ripley=F)
.RA
.AG object
result of a model fit using the `survreg' function.
.OA
.AG newdata
data for prediction.  If absent predictions are for the
subjects used in the original fit.
.AG type
the type of predicted value. 
This can be on the original scale of the data (response),
the linear predictor (`"linear"', with `"lp"' as an allowed abbreviation),
a predicted quantile on the original scale of the data (quantile),
a quantile on the linear predictor scale (uquantile),
or the matrix of terms for the linear predictor (terms).
At this time `"link"' and linear predictor (`"lp"') are identical.
.AG se.fit
if TRUE, include the standard errors of the prediction in the result.
.AG terms
subset of terms.  The default for residual type`terms' is a matrix with
one column for every term (excluding the intercept) in the model.
.AG p
vector of percentiles.  This is used only for quantile predictions.
.AG ripley
temporary arg. Compute the standard errors of quantile predictions in
the way shown in an a draft of Ripley and Venables, i.e., partially ignoring
the variance in the estimate of scale.
.RT
a vector or matrix of predicted values.
.SH REFERENCES
Escobar and Meeker (1992). Assessing influence in regression analysis with
censored data. \fIBiometrics,\fP 48, 507-528.
.SA
`survreg', `residuals.survreg'
.EX
# Draw figure 1 from Escobar and Meeker
fit <- survreg(Surv(time,status) ~ age + age^2, data=stanford2,
	dist='lognormal')
plot(stanford2$age, stanford2$time, xlab='Age', ylab='Days',
	xlim=c(0,65), ylim=c(.01, 10^6), log='y')
pred <- predict(fit, newdata=list(age=1:65), type='quantile',
	         p=c(.1, .5, .9))
matlines(1:65, pred, lty=c(2,1,2), col=1)
.KW survival
.WR
