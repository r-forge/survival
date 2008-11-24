#
# This was to test the new rate tables, after I had changed from mdy.date()
#  to julian().  This test doesn't have to be repeated again.
# I know that position "data" contains the rate tables from StatSci.
# 
# Later note: the Statsci rate tables are way out of date.  They ignored the
#  1996 update of the data!  So I'll have to compare to Mayo Splus3.3 version.
#
# The test data purposefully covers a large range of years
tdata <- data.frame(age=1:100 * 365, 
                    year=julian(sample(1:12, 100, replace=T),
                                rep(15,100),
                                sample(1930:2000, 100, replace=T)),
                    sex = rep(1:2, 50),
                    race= rep(1:2, length=100))

tablelist <- c('survexp.az', 'survexp.azr', 'survexp.fl', 'survexp.flr',
               'survexp.mn', 'survexp.mnwhite', 'survexp.us',
               'survexp.usr', 'survexp.wnc')


# This part has to be run in Splus3.4
ntable <- length(tablelist)
rmat <- matrix(0, nrow=50, ncol=ntable)
for (i in 1:ntable) {
    fit <- survexp(~1, data=tdata, ratetable=get(tablelist[i]),
                    times= 200*(1:50))
    rmat[,i] <- fit$surv
    }
data.dump('rmat')

#
# This part in 6.0
#  (remember to copy tdata over, not re-make it!)
#

for (i in 1:ntable) {
    fit2 <- survexp(~1, data=tdata, ratetable= get(tablelist[i]),
                    times= 200*(1:50))
    xx <- all.equal(rmat[,i], as.vector(fit2$surv))
    if (!is.logical(xx) || !xx)
        cat("Problem with ", tablelist[i], "\n")
    else cat(i, " ok\n")       
    }

rm(tdata, tablelist, fit1, fit2, xx)

