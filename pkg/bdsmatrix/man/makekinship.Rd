\name{makekinship}
\alias{makekinship}
\title{Create a sparse kinship matrix}
\description{
Compute the overall kinship matrix for a collection of families, and
store it efficiently.
}
\usage{makekinship(famid, id, father.id, mother.id, unrelated=0)}
\arguments{
  \item{famid}{
    a vector of family identifiers
  }
  \item{id}{
    a vector of unique subject identifiers
  }
  \item{father.id }{
    for each subject, the identifier of their biolgical father
  }
  \item{mother.id }{
    for each subject, the identifier of thier biological mother
  }
  \item{unrelated}{
    subjects with this family id are considered to be unrelated singletons,
    i.e., not related to each other or to anyone else.
  }
}
\value{
  a block-diagonal-symmetric kinship matrix of class \code{bdsmatrix}
}
\details{
  For each family of more than one member, the 
  \code{kinship} function is called to calculate a
  per-family kinship matrix.
  These are stored in an efficient way into a single block-diagaonal sparse
  matrix object, taking advantage of the fact that between family entries
  in the full matrix are all 0.
  Unrelated individuals are considered to be families of size 0, and are
  placed first in the matrix.

  The final order of the rows within this matrix will not necessarily be
  the same as in the origianl data, since each family must be contiguous.
  The dimnames of the matrix contain the id variable for each row/column.
  Also note that to create the kinship matrix for a subset of the data
  it is necessary to create the full kinship
  matrix first and then subset it.  
  One cannot first subset the data and then call the function.
  For instance, a call using only the female data would not detect that a
  particular man's sister and his daughter are related.
}
\seealso{\code{\link{kinship}}, \code{\link{makefamid}}}
\examples{
# Data set from a large family study of breast cancer
#  there are 26050 subjects in the file, from 426 families
table(mnbreast$sex)
length(unique(mnbreast$famid))
kin1 <- makekinship(mnbreast$famid, mnbreast$id,
                    mnbreast$dadid, mnbreast$momid)
dim(kin1)
class(kin1)

# The next line shows that less than 2\% of the elements are non-zero
length(kin1@blocks)/ prod(dim(kin1))

# kinship matrix for the females only
femid <- mnbreast$id[mnbreast$sex=='F']
femindex <- !is.na(match(dimnames(kin1)[[1]], femid))
kin2 <- kin1[femindex, femindex]
#
# Note that "femindex <- match(femid, dimnames(kin1)[[1]])" is wrong, since
#  then kin1[femindex, femindex] might improperly reorder the rows/cols 
# However sort(match(femid, dimnames(kin1)[[1]])) would be okay.
}
\keyword{genetics}

