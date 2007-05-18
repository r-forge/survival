\name{rowsum}
\alias{rowsum}
\title{
  Row Sums of a Matrix, Based on a Grouping Variable. 
}
\description{
Compute sums across rows of a matrix for each level of a grouping variable. 
}
\usage{
rowsum(x, group)
}
\arguments{
\item{x}{
 a matrix or vector of numeric data.  Missing values are allowed. 
}
\item{group}{
 a vector giving the grouping, with one element per row of \code{x}. 
Missing values are not allowed. 
}
}
\value{
a matrix containing the sums.  There will be one row per unique value 
of \code{group} and the same number of columns as \code{x}. 
}
\seealso{
\code{\link{tapply}}. 
}
\examples{
x <- matrix(runif(100), ncol=5) 
group <- sample(1:8, 20, T) 
xsum <- rowsum(x, group) 
# same result another way, slower, and temp may be much larger than x 
temp <- model.matrix( ~ a - 1, data.frame(a=as.factor(group))) 
xsum2 <- t(temp) \%*\% x 
}
\keyword{manip}
\keyword{iteration}
\keyword{category}
% docclass is function
% Converted by Sd2Rd version 37351.
