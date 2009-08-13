# Automatically generated from all.nw using noweb
expand.nested <- function(x) {
    x[[1]] <- as.factor(x[[1]])[,drop=T]
    if (length(x) >1) {
        for (i in seq(2, length(x), by=1)) {
            x[[i]] <- strata(x[[i-1]], x[[i]], shortlabel=TRUE, sep='/')
            }
       } 
    x
    }
