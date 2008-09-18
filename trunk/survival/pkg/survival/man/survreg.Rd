\name{survreg}
\alias{survreg}
\title{
  Regression for a Parametric Survival Model 
}
\description{
Fit a parametric survival regression model. 
The function \code{censorReg} is an alternative function for fitting these type 
of model. 
}
\usage{
survreg(formula=formula(data), data=sys.parent(), subset, 
        na.action, link="log", dist="extreme", init=NULL, 
        fixed=list(), control, model=F, x=F, y=T, ...)
}
\arguments{
\item{formula}{
a formula expression as for other regression models. 
The response is usually a survival object as returned by the \code{Surv} function. 
See the documentation for \code{Surv}, \code{lm} and \code{formula} for details. 
}
\item{data}{
a data frame in which to interpret the variables named in 
the \code{formula}, or the \code{subset} argument. 
}
\item{subset}{
expression saying that only a subset of the rows of the data 
should be used in the fit. 
}
\item{na.action}{
a missing-data filter function, applied to the model.frame, after any 
\code{subset} argument has been used.  Default is \code{options()\$na.action}. 
}
\item{link}{
a character string specifying the transformation to be used on the y variable. 
Possible values are \code{"log"} and \code{"identity"}. 
}
\item{dist}{
a character string specifying the assumed distribution for the transformed y 
variable. 
Possible values are \code{"extreme"}, \code{"logistic"}, \code{"gaussian"}, \code{"exponential"}, 
and \code{"rayleigh"}. 
}
\item{init}{
a vector of initial values for the parameters. 
}
\item{fixed}{
a list of fixed parameters, most often just the scale. 
}
\item{control}{
a list of three parameters (\code{maxiter}, \code{rel.tolerance} and \code{failure}) 
returned by \code{survreg.control}. 
See \code{survreg.control} and the EXAMPLES section for details. 
}
\item{model}{
a logical value, if \code{TRUE} then the model frame is included as component 
\code{model} in the object returned by the function. 
}
\item{x}{
a logical value, if \code{TRUE} then the model matrix is returned as component 
\code{x} in the object returned by the function. 
}
\item{y}{
a logical value, if \code{TRUE} then the response is returned as component \code{y} 
in the object returned by the function. 
}
\item{\dots}{
all the optional arguments to \code{lm}, including \code{singular.ok}. 
}
}
\value{
an object of class \code{survreg} is returned, which inherits from class \code{glm}. 
See \code{survreg.object} for details. 
}
\section{Details}{
This routine is not as robust against nearly singular X matrices as 
\code{lm}; the problem occurs when the covariance matrix is explicitly inverted 
using the \code{solve} function.  This can sometimes be solved by subtracting 
the mean from all continuous covariates. 
}
\seealso{
\code{\link{censorReg}},  \code{\link{formula}},  \code{\link{lm}},  \code{\link{Surv}},  \code{\link{survreg.control}},  \code{\link{survreg.object}}.   
}
\examples{
# Fit an exponential model 
survreg(Surv(futime, fustat) ~ ecog.ps + rx, data=ovarian, 
        dist='extreme', link='log', fixed=list(scale=1))   
# Using survreg.control to supply the control argument 
survreg(Surv(days, event) ~ voltage, data=capacitor, 
        control=survreg.control(maxiter = 40)) 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
