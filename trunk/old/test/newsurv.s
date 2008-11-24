#
# Some checks on the new survival tables
#
all.equal(survexp.oldus[-1,,1:3], survexp.us[-(1:4),,3:5])

# First year total hazard
temp1 <- c(1,6,21, 365.24-28) %*% survexp.us[1:4,1,3:5]
temp2 <- 365.24*survexp.oldus[1,1,1:3]
all.equal(as.vector(temp1), as.vector(temp2)) 

#
# US by race
#
temp1 <- aperm(survexp.usr[5:113,,,3:5], c(1,2,4,3))
temp2 <- survexp.oldusr[2:110,,1:3,]
all.equal(as.vector(temp1), as.vector(temp2))

# First year total hazard
for (i in 1:3) {
    temp1 <- c(1,6,21, 365.24-28) %*% survexp.usr[1:4,1,i, 3:5]
    temp2 <- 365.24*survexp.oldusr[1,1,1:3,i]
    print(all.equal(as.vector(temp1), as.vector(temp2))) 
    }

# Black 1950-60 is a duplicate
all.equal(as.vector(survexp.usr[,,'nonwhite', c("1950", "1960")]),
	  as.vector(survexp.usr[,,'black',    c("1950", "1960")]))

# West North Central
all.equal(survexp.wnc[,,1:8], survexp.oldwnc[,,1:8])
all.equal(survexp.wnc[-(1:2),,7:10], survexp.mnwhite[-1,,3:6])

# Minnesota
all.equal(survexp.mn[,,1:2], survexp.oldmn[,,1:2])

# Florida
#  I didn't save old ones, and the Statsci ones differ on not having
#   a class "date" for the year attributes.  So use as.vector to compare
temp <- get('survexp.fl', where='data')
all.equal(as.vector(survexp.fl[,,1:2]), as.vector(temp[,,1:2]))

