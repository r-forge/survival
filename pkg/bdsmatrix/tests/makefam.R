library(bdsmatrix)
#
# A pedigree that caused problems for an earlier version of the
#  code due to a double marriage.  It mistakenly saw 2 famiies.
id <- 1:20
mom<- c(0,0,0,2,2,2,0,2,0, 0,2,2,0,2,0,2, 7,7, 11,14)
dad<- c(0,0,0,1,1,1,0,1,0, 0,3,3,0,3,0,3, 8,8, 10,13)
sex<- c(0,1,0,0,1,1,1,0,1, 0,1,0,0,1,1,1, 0,0, 1, 0)

temp<- makefamid(id, dad, mom)
temp2 <- rep(1,20)
temp2[c(9, 15)] <- 0  # the only 2 marry-ins
all.equal(temp, temp2)

