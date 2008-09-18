\name{summary.survfit}
\alias{summary.survfit}
\title{
  Summary of a Survival Curve 
}
\description{
Returns a list containing the survival curve, confidence limits for the 
curve, and other information. 
}
\usage{
summary.survfit(fit, times=&lt;&lt;see below&gt;&gt;, censored=F, scale=1, extend=F)
}
\arguments{
\item{fit}{
the result of a call to the \code{survfit} function. 
}
\item{times}{
vector of times; 
the returned matrix will contain 1 row for each time. 
This must be in increasing order and missing values are not allowed. 
If \code{censored=T}, the default \code{times} vector contains all the unique times in 
\code{fit}, 
otherwise the default \code{times} vector uses only the event (death) times. 
}
\item{censored}{
logical value: if TRUE, the censoring times are included in the output. 
This is ignored if the \code{times} argument is present. 
}
\item{scale}{
numeric value to rescale the survival time, e.g., if the input data to 
\code{survfit} were in 
days, \code{scale = 365.25} would scale the output to years. 
}
\item{extend}{
logical value: if TRUE, prints information for all specified \code{times}, 
even if there are no subjects left at the end of the specified \code{times}. 
This is only valid if the \code{times} 
argument is present. 
}
}
\value{
a list with the following components: 

\item{surv}{
the estimate of survival at time t+0. 
}
\item{time}{
the timepoints at which the curve has a step. 
}
\item{n.risk}{
the number of subjects at risk at time t-0 
(but see the comments on weights in the \code{survfit} help file). 
}
\item{n.event}{
if the \code{times} argument is missing, this column is the number of 
events that occurred at time t. 
Otherwise, it is the cumulative number of events that have occurred 
since the last time listed until time t+0. 
}
\item{n.entered}{
if the \code{type} argument is \code{counting} and the \code{times} argument is 
missing, this column is the number of subjects that entered at time t. 
Otherwise, it is the cumulative number of subjects that have entered 
since the last time listed until time t+0.  This is only valid if the 
\code{type} argument is \code{counting}. 
}
\item{n.exit.censored}{
if the \code{type} argument is \code{counting} and the \code{times} argument is 
missing, this column is the number of subjects that left without an 
event at time t. 
Otherwise, it is the cumulative number of subjects that have left 
without an event 
since the last time listed until time t+0.  This is only valid if the 
\code{type} argument is \code{counting}. 
}
\item{std.err}{
the standard error of the survival value. 
}
\item{conf.int}{
level of confidence for the confidence intervals of survival. 
}
\item{lower CI}{
lower confidence limits for the curve. 
}
\item{upper CI}{
upper confidence limits for the curve. 
}
\item{strata}{
indicates stratification of curve estimation.  If \code{strata} is not \code{NULL}, 
there are multiple curves in the result and the \code{surv}, \code{time}, \code{n.risk}, etc.  
vectors will contain multiple curves, pasted end to end.  
The levels of \code{strata} (a factor) are the labels for the curves. 
}
\item{call}{
the statement used to create the \code{fit} object. 
}
\item{na.action}{
same as for \code{fit}, if present. 
}
\item{table}{
table of information that is returned from \code{print.survfit} function. 
}
\item{type}{
type of data censoring.  Passed from \code{survfit} function. 
}
}
\seealso{
\code{\link{survfit}},  \code{\link{print.summary.survfit}}.   
}
\examples{
summary(survfit(Surv(time, status) ~ group, data = leukemia))  
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
