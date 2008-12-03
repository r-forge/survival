#
# Generate each of the messages from is.ratetable
#
temp <- runif(21*2*4)

# Good
attributes(temp) <- list(dim=c(21,2,4),
    dimnames=list(c(as.character(75:95)), c("male","female"),
                  c(as.character(2000:2003))),
    dimid=c("age","sex","year"),
    factor=c(0,1,0),
    cutpoints=list(c(75:95), NULL, chron('1/1/2000') +c(0:3)*366.25),
    class='ratetable')
is.ratetable(temp)

# Factor problem + cutpoints length
attributes(temp) <- list(dim=c(21,2,4),
    dimnames=list(c(as.character(75:95)), c("male","female"),
                  c(as.character(2000:2003))),
    dimid=c("age","sex","year"),
    factor=c(1,1,0),
    cutpoints=list(c(75:95), NULL, chron('1/1/2000') +c(0:4)*366.25),
    class='ratetable')
is.ratetable(temp, verbose=T)
 
                    
# missing dimid attribute + unsorted cutpoint
attributes(temp) <- list(dim=c(21,2,4),
    dimnames=list(c(as.character(75:95)), c("male","female"),
                  c(as.character(2000:2003))),
    factor=c(0,1,0),
    cutpoints=list(c(75:95), NULL, chron('1/1/2000') +c(4:1)*366.25),
    class='ratetable')
is.ratetable(temp, verbose=T)

# wrong dimname and dimid
attributes(temp) <- list(dim=c(21,2,4),
    dimnames=list(c(as.character(75:95)), c("male","female"),
                  c(as.character(2000:2003))),
    dimid=c("age","sex","year", "zed"),
    factor=c(0,1,0),
    cutpoints=list(c(75:95), NULL, chron('1/1/2000') +c(0:3)*366.25),
    class='ratetable')
is.ratetable(temp, verbose=T)
