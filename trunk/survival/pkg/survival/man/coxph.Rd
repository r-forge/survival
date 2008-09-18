\name{coxph}
\alias{coxph}
\title{
  Fit Proportional Hazards Regression Model 
}
\description{
Fits a Cox proportional hazards regression model. 
Time dependent variables, time dependent strata, multiple events per subject, 
and other extensions are incorporated using the counting process formulation 
of Andersen and Gill. 
}
\usage{
coxph(formula, data=sys.parent(), weights, subset, 
      na.action, eps = 1e-07, init, iter.max=10,  
      method=c("efron","breslow","exact"), 
      singular.ok=T, robust=F, 
      model=F, x=F, y=T)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a `~' operator, and 
the terms on the right.  The response must be a survival object as 
returned by the \code{Surv} function. 
}
\item{data}{
a \code{data.frame} in which to interpret the variables named in 
the \code{formula}, or in the \code{subset} and the \code{weights} argument. 
}
\item{weights}{
vector of case weights.  If \code{weights} is a vector of integers, then the 
estimated coefficients are equivalent to estimating the model from data 
with the individual cases replicated as many times as indicated by 
\code{weights}. 
}
\item{subset}{
expression indicating which subset of the rows of data should be used in 
the fit.  This can be a logical vector (which is replicated to have 
length equal to the number of observations), a numeric vector indicating 
which observation numbers are to be included (or excluded if negative), 
or a character vector of row names to be included.  All observations are 
included by default. 
}
\item{na.action}{
a missing-data filter function.  This is applied to the \code{model.frame} after any 
subset argument has been used.  Default is \code{options()\$na.action}. 
}
\item{eps}{
convergence threshold.  Iteration will continue until the relative change 
in the log-likelihood is less than \code{eps}.  Default is 1e-07. 
}
\item{init}{
vector of initial values of the iteration.  Default initial 
value is zero for all variables. 
}
\item{iter.max}{
maximum number of iterations to perform.  Default is 10. 
}
\item{method}{
a character string specifying the method for tie handling.  If there  
are no tied death times all the methods are equivalent. 
Nearly all Cox regression programs use the Breslow method by default, 
but 
not this one. 
The Efron approximation is used as the default here, as it is much more 
accurate when dealing with tied death times, and is as efficient 
computationally. 
The exact method computes the exact partial likelihood, which is 
equivalent 
to a conditional logistic model.  If there are a large number of ties 
the 
computational time will be excessive. 
}
\item{singular.ok}{
logical value indicating how to handle collinearity in the model matrix. 
If \code{TRUE}, the program will automatically skip over columns of the X 
matrix that are linear combinations of earlier columns.  In this case the 
coefficients for such columns will be NA, and the variance matrix will 
contain zeros. For ancillary calculations, such as the linear predictor, the 
missing coefficients are treated as zeros. 
}
\item{robust}{
if \code{TRUE} a robust variance estimate is returned.  Default is \code{TRUE} if the 
model includes a \code{cluster} operative, \code{FALSE} otherwise. 
}
\item{model}{
logical value: if \code{TRUE}, the model frame is returned in component model 
}
\item{x}{
logical value: if \code{TRUE}, the model frame is returned in component \code{x}. 
}
\item{y}{
logical value: if \code{TRUE}, the model frame is returned in component \code{y}. 
}
}
\value{
an object of class \code{"coxph"} representing the fit. 
See \code{coxph.object} for details.  
}
\section{Side Effects}{
Depending on the call, the \code{predict}, \code{residuals}, and \code{survfit} routines may 
need to reconstruct the x matrix created by \code{coxph}.  Differences in the 
environment, such as which data frames are attached or the value of 
\code{options()\$contrasts}, may cause this computation to fail or worse, to be 
incorrect.  See the survival overview document for details. 
}
\details{
The proportional hazards model is usually expressed in terms of a 
single survival time value for each person, with possible censoring. 
Andersen and Gill reformulated the same problem as a counting process; 
as time marches onward we observe the events for a subject, rather 
like watching a Geiger counter. 
The data for a subject is presented as multiple rows or "observations", 
each 
of which applies to an interval of observation (start, stop]. 
}
\section{Special terms}{
There are two special terms that may be used in the model equation. 
A \code{strata} term identifies a stratified Cox model; separate baseline 
hazard 
functions are fit for each strata. 
The \code{cluster} term is used to compute a robust variance for the model. 
The term \code{+ cluster(id)}, where \code{id == unique(id)}, is equivalent to 
specifying the \code{robust=T} argument, and produces an approximate 
jackknife estimate of the variance.  If the \code{id} variable were not 
unique, but instead 
identifies clusters of correlated observations, then the variance 
estimate 
is based on a grouped jackknife. 
}
\section{Convergence}{
In certain data cases the actual MLE estimate of a 
coefficient is infinity, e.g., a dichotomous variable where one of the 
groups has no events.  When this happens the associated coefficient 
grows at a steady pace and a race condition will exist in the fitting 
routine: either the log likelihood converges, the information matrix 
becomes effectively singular, an argument to exp becomes too large for 
the computer hardware, or the maximum number of interactions is exceeded. 
The routine attempts to detect when this has happened, not always successfully. 
}
\section{References}{
Andersen, P. and Gill, R. (1982). 
Cox's regression model for counting processes, a large sample study. 
\emph{Annals of Statistics}
\bold{10}, 1100-1120. 

Therneau, T., Grambsch, P., and Fleming. T. (1990). 
Martingale based residuals for survival models. 
\emph{Biometrika}
\bold{77}, 147-160. 
}
\seealso{
\code{\link{cluster}},  \code{\link{strata}},  \code{\link{Surv}},  \code{\link{survfit}}.   
}
\examples{
# Create the simplest test data set 
test1 <- list(time=c(4,3,1,1,2,2,3), 
              status=c(1,1,1,0,1,1,0), 
              x=c(0,2,1,1,1,0,0), 
              sex=c(0,0,0,0,1,1,1)) 
# Fit a stratified model 
coxph(Surv(time, status) ~ x + strata(sex), test1) 
# Create a simple data set for a time-dependent model 
test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 
summary(coxph(Surv(start, stop, event) ~ x, test2)) 
# Fit a stratified model, clustered on patients 
bladder1 <- bladder 
bladder1$start <- NULL 
bladder1 <- bladder1[bladder1$enum < 5, ] 
coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
      cluster(id), bladder1, method='breslow') 
}
\keyword{models}
\keyword{regression}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
