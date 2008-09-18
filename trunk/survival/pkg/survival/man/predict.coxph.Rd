\name{predict.coxph}
\alias{predict.coxph}
\title{
  Make predictions from a fitted coxph object. 
}
\description{
Creates predicted values from a coxph fit. 
}
\usage{
predict.coxph(object, newdata, type=c("lp", "risk", "expected", "terms"),  
se.fit=F, terms=labels.lm(object), collapse)
}
\arguments{
\item{object}{
the results of a coxph fit. 
}
\item{newdata}{
data for prediction.  If absent predictions are for the 
subjects used in the original fit. 
}
\item{type}{
the type of predicted value. 
Choices are the linear predictor (lp), the risk score exp(lp) (risk), 
the expected number of events given the covariates and follow-up time 
(expected), and the terms of the linear predictor (terms). 
}
\item{se.fit}{
if TRUE, pointwise standard errors are produced for the predictions. 
}
\item{terms}{
if type="terms", this argument can be used to specify which terms should be 
included; the default is all. 
}
\item{collapse}{
optional vector of subject identifiers.  
If specified, the output will contain one entry per subject rather than one 
entry per observation. 
}
}
\value{
a vector or matrix of predictions, or a list containing the predictions 
(element "fit") and their standard errors (element "se.fit") if the se.fit 
option is TRUE. 
}
\details{
This function is a method for the generic function \code{predict} for class 
\code{coxph}.  It can be invoked by calling \code{predict} for an object of the 
appropriate class, or by calling predict.coxph directly. 
}
\seealso{
coxph, residuals.coxph  
}
\examples{
> options(na.action='na.omit') 
> fit <- coxph(Surv(time, status) ~ age + ph.ecog + strata(inst), lung) 
> mresid <- lung$status - predict(fit, type='expected') #Martingale resid 
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
