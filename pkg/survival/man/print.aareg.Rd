\name{print.aareg}
\alias{print.aareg}
\title{
  Print an aareg object
}
\description{
Print out a fit of Aalen's additive regression model
}
\usage{
\method{print}{aareg}(x, maxtime, test=c("aalen", "nrisk"))
}
\arguments{
\item{x}{
the result of a call to the \code{aareg} function
}
\item{maxtime}{
the upper time point to be used in the test for non-zero slope
}
\item{test}{
the weighting to be used in the test for non-zero slope.
The default weights are based on the variance of each coefficient, as
a function of time.  The alternative weight is proportional to the number
of subjects still at risk at each time point.
}
}
\value{
the calling argument is returned.
}
\section{Side Effects}{
the results of the fit are displayed.
}
\details{
The estimated increments in the coefficient estimates can become quite
unstable near the end of follow-up, due to the small number of observations
still at risk in a data set.
Thus, the test for slope will sometimes be more powerful if this last
`tail' is excluded.
}
\section{References}{
Aalen, O.O. (1989). A linear regression model for the analysis of life times.
Statistics in Medicine, 8:907-925.
}
\seealso{
aareg
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
