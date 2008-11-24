#
# Set up for the test
#
#options(echo=T)
options(na.action="na.exclude", 
	contrasts=c('contr.treatment','contr.poly'), conflicts.ok=T)
attach("..")
date()

aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
