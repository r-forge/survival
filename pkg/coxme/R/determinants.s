aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# This is an exercise to help me understand determinants
#   Start with a symmetric positive definite matrix, sparse first column
#
temp <- matrix(c(5,1,0,0,1,6,2,2,0,2,7,3,0,2,3,8), ncol=4)

# get the chol and the svd
#   The determinant should be the product of the diagonal of the 
#  general chol, and of the d matrix from svd
gtemp <- gchol(temp)
stemp  <- svd(temp)
aeq(prod(diag(gtemp)), prod(stemp$d))

# For the upper 3x3 columns, the original gchol still works
stemp3 <- svd(temp[1:3,1:3])
aeq(prod(diag(gtemp)[1:3]), prod(stemp3$d))

#
# Ok, now what is the determinant of the upper 3x3 of the inverse?
#
