\name{print.survfit}
\alias{print.survfit}
\title{
  Print a Short Summary of a Survival Curve 
}
\description{
Print number of observations, number of events, the mean survival and its  
standard error, and the median survival with confidence limits for the  
median. 
}
\usage{
print.survfit(x, scale=1)
}
\arguments{
\item{x}{
the result of a call to the \code{survfit} function. 
}
\item{scale}{
a numeric value to rescale the survival time, e.g., 
if the input data to survfit were in days, 
\code{scale=365} would scale the printout to years. 
}
}
\value{
x, with the invisible flag set to prevent printing. 
}
\section{Side Effects}{
the number of observations, the number of events, the mean survival and its 
standard error, and the median survival with its confidence interval are 
printed.  If there are multiple curves, there is one line of output for each. 
}
\details{
The mean and its variance are based on a truncated estimator.  That is, if the 
last observation(s) is not a death, then the survival curve estimate does not 
go to zero and the mean is undefined.  In such a case, the estimator is based 
on an assumption that the true curve goes to zero just beyond the last 
observed follow up time; it will systematically underestimate the true mean. 

The median and its confidence interval are defined by drawing a horizontal 
line at 0.5 on the plot of the survival curve and its confidence bands. 
The intersection of the line with the lower CI band defines the lower limit 
for the median's interval, and similarly for the upper band.  If any of the 
intersections is not a point, then we use the smallest point of intersection, 
e.g., if the survival curve were exactly equal to 0.5 over an interval. 
}
\section{References}{
Miller, Rupert G., Jr. (1981).  
\emph{Survival Analysis.}
New York:Wiley, p 71. 
}
\seealso{
\code{\link{summary.survfit}}. 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
