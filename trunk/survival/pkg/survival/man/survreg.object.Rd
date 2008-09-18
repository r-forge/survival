\name{survreg.object}
\alias{survreg.object}
\title{
  Parametric Survival Model Object 
}
\description{
Classes of objects that result from fitting a parametric survival model. 
}
\section{Generation}{
This class of objects is returned by the \code{survreg} function 
to represent a fitted parametric survival model. 
Class \code{survreg} inherits from class \code{glm}, since it is fit by iterative 
reweighted least squares; the object returned has all the components of a 
weighted least squares object. 
}
\section{Methods}{
Objects of this class have methods for the functions \code{print}, 
\code{summary}, \code{predict}, and \code{residuals}. 
}
\section{Structure}{
The following components must be included in a legitimate \code{survreg} object. 
The residuals, fitted values, coefficients and effects should be extracted 
by the generic functions of the same name, rather than 
by the \code{"\$"} operator.  
}
\arguments{
\item{coefficients}{
the coefficients of the \code{linear.predictors}, which multiply  the 
columns of the model matrix.  It does not include the estimate of error  
(sigma).  The names of the coefficients are the names of the 
single-degree-of-freedom effects (the columns of the model matrix). 
If the model is over-determined there will be missing values in the  
coefficients corresponding to nonestimable coefficients. 
}
\item{family}{
a 2 element character vector giving the name of the family and 
the link; mainly for printing purposes. 
}
\item{parms}{
the parameters of the model that are not coefficients of the model matrix. 
The first of these will always be \code{log(sigma)}. 
}
\item{fixed}{
a vector the same length as \code{parms}, where 1 indicates a parameter that 
was fixed at its starting value and was not part of the iteration. 
}
\item{var}{
the final variance matrix, including both coefficients and free parameters. 
}
\item{dresiduals}{
the deviance residuals. 
}
\item{deviance}{
minus twice the difference between the maximized log-likelihood under the 
fitted model and a saturated model. Similar to the residual sum of squares. 
}
\item{null.deviance}{
the deviance corresponding to the model with only an intercept term, and 
with \code{parms} fixed at their final values. 
}
\item{loglik}{
the log-likelihood for the final model. 

The object will also have the components of a \code{glm} object not mentioned 
above: \code{linear predictors}, \code{fitted.values}, \code{residuals}, \code{effects}, \code{R},  
\code{rank}, \code{assign}, \code{df.residual}, \code{weights}, \code{call}, \code{iter}, \code{contrasts},  
\code{terms} and \code{formula}.  See \code{glm.object}. 
}
}
\seealso{
\code{\link{glm.object}},  \code{\link{lm.object}},  \code{\link{survreg}}.   
}
\keyword{regression}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
