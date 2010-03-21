\name{kindepth}
\Rdversion{1.1}
\alias{kindepth}
\title{Compute pedigree depth}
\description{For each subject in a pedigree, compute their depth.
  This is the number of generations of ancestors they have in the
  pedigree.}
\usage{kindepth(id, dad.id, mom.id, align = FALSE)}
\arguments{
  \item{id}{Subject identifier}
  \item{dad.id}{Identifier of the biological father}
  \item{mom.id}{Identifier of the biological mother}
  \item{align}{If align=T, go one step further and try to make both
    parents of each child have the same depth.  (This is not always
    possible). 
    }
}
\details{
  Depth is defined to be 0 for a founder, and 1+ max(father's depth,
  mother's depth) otherwise.
  The basic algorithm is to make founder's depth=0, children of founders
  =1, children of children of founders =2, etc.  
  
  With the \code{align} argument, some founders may be given values
  greater than zero, the final value is each subject's distance from the
  top of the page in a `nice' drawing of the pedigree.
  (It mostly affects pedigress that "join in the middle" via a marriage.)}

\value{
  an integer vector giving the depth of each subject.}
\author{Terry Therneau}
\seealso{\code{\link{kinship}}}
\keyword{ genetics }
