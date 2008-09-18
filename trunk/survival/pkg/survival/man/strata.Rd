\name{strata}
\alias{strata}
\title{
  Identify Strata Variables 
}
\description{
This is a special function used in the context of the Cox survival model. 
It identifies stratification variables when they appear on the right hand 
side of a formula. 
}
\usage{
strata(..., na.group=F, sep=', ')
}
\arguments{
\item{\dots}{
any number of variables.  All must be the same length. 
}
\item{na.group}{
a logical variable, if TRUE, then missing values are treated as a 
distinct level of each variable. 
}
\item{sep}{
the character used to separate groups, in the created label
}
}
\value{
a new factor, whose levels are all possible combinations of the factors 
supplied as arguments. 
}
\details{
The result is identical to the \code{interaction} function, 
but for the labeling of the factors (\code{strata} is more verbose). 
}
\seealso{
\code{\link{coxph}}. 
}
\examples{
coxph(Surv(futime, fustat) ~ age + strata(rx), data=ovarian) 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
