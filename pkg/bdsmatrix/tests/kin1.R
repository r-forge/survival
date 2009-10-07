library(bdsmatrix)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# The test case from Lange, chapter 5
#
id <- 1:6
momid <- c(0,0,2,2,4,4)
dadid <- c(0,0,1,1,3,3)

xx <- kinship(id, dadid, momid)
aeq(xx, c(4,0,2,2,2,2, 0,4,2,2,2,2, 2,2,4,2,3,3, 
          2,2,2,4,3,3, 2,2,3,3,5,3, 2,2,3,3,3,5) /8)


# And here is an an odd one with cross marriages, but no inbreeding
#
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))

xx <- kinship(test1$id, test1$dad, test1$mom)
all.equal(xx, t(xx))
aeq(diag(xx), rep(.5,14))

aeq(8*xx[1,], c(4,0,0,0,2,2,0,0, 1,0,0,0, 0,0))
aeq(8*xx[2, 3:14], c(0,0,2,2,0,0, 1,2,0,0,0, 1))
aeq(8*xx[3, 4:14], c(0,0,0,2,2, 2,1,0,0,0,.5))
aeq(8*xx[4, 5:14], c(0,0,2,2, 0,1, 0,0,0, .5))
aeq(8*xx[5, 6:14], c(2,0,0,1, 1,0,0,0,.5))
aeq(8*xx[6, 7:14], c(0,0,2,1, 0,0,0, .5))
aeq(8*xx[7, 8:14], c(2,1, 2, 0,0,0, 1))
aeq(8*xx[8, 9:14], c(1,1, 0,0,0,.5))
aeq(8*xx[9, 10:14], c(1, 0,0,0, .5))
aeq(8*xx[10,11:14], c(0,0,0,2))
aeq(8*xx[11, 12:14], c(0,2,1))
aeq(8*xx[12,13:14], c(2,1))
aeq(xx[13,14], .25)

