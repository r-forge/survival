# $Id: summary.ratetable.s,v 1.3 2006-08-28 15:40:30 m015733 Exp $
#
# Print out information about a rate table: it's dimensions and keywords
#
summary.ratetable <- function(x, ...) {  
  alist <- attributes(x)
  ndim <- length(alist$dim)
  if(ndim > 1)
    cat("Rate table with dimensions:\n")
  else cat("Rate table with 1 dimension:\n")
  fac <- alist$factor
  temp1 <- format(alist$dimid)
  for(i in 1:ndim) {
    cat("    ", temp1[i], ": ", sep = "")
    if(fac[i] == 1) {
      cat("discrete factor with legal values of (",
          paste(alist$dimnames[[i]], collapse = ", "), ")\n", sep = "")
    }
    else if(fac[i] == 0)
      cat("time variable with", format(alist$dim[i]), 
          "categories\n")
    else cat("time variable with", format(alist$dim[i]), 
             "categories (interpolated)\n")
  }
}
