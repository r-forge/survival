\name{makefamid}
\alias{makefamid}
\title{Identify family groups}
\description{
Given a set of parentage relationships, this subdivides a set of
subjects into families. 
}
\usage{makefamid(id, father.id, mother.id)}
\arguments{
  \item{id}{a vector of unique subject identifiers}
  \item{father.id }{for each subject, the identifier of their biolgical father}
  \item{mother.id }{for each subject, the identifier of thier biological mother}
}

\value{
  a vector of family identifiers.
  Individuals who are not blood relatives of anyone else in the data
  set are assigned a family id of 0.
}
\details{
  This function may be useful to create a family identifier if none
  exists in the data (rare),
  to check for anomalies in a given family identifier
  (see the \code{familycheck} function),
  or to create a more space and time efficient kinship matrix by
  separating out marry-ins without children as 'unrelated'.
}
\seealso{
  makefamid, kinship, makekinship
}
\examples{
newid <- makefamid(mnbreast$id, mnbreast$dadid, mnbreast$momid)
table(newid==0)  # Nearly 1/3 are not blood relatives (newid==0)

kin1 <- makekinship(mnbreast$famid, mnbreast$id, mnbreast$dadid,
                                    mnbreast$momid)
kin2 <- makekinship(newid, mnbreast$id, mnbreast$dadid,
                    mnbreast$momid, unrelated=0)
dim(kin2)  #both will be 26050 square sparse matrices
dim(kin1)

length(kin2@blocks)/length(kin1@blocks)  #Ratio of non-zero elements
# Basing kin1 on newid rather than mnbreast$famid (where marry-ins were each
#   labeled as members of one of the 426 families) reduced its size by just
#   less than half.
}
\keyword{genetics}
