\name{cluster}
\alias{cluster}
\title{
  Identify clusters. 
}
\description{
This is a special function used in the context of the Cox model.  It 
identifies correlated groups of observations, and is used on the right hand 
side of a formula. 
}
\usage{
cluster(x)
}
\arguments{
\item{x}{
A character, factor, or numeric variable. 
}
}
\value{
x 
}
\details{
The function's only action is semantic, to mark a variable as the 
cluster indicator. 
}
\seealso{
\code{\link{coxph}},  \code{\link{Surv}}.   
}
\examples{
coxph(Surv(futime, fustat) ~ age + cluster(rx), data = ovarian) 
}
\keyword{survival4}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 37351.
