## $Id: print.ratetable.s,v 5.1 2006-02-23 23:03:26 lunde Exp $

print.ratetable <- function(x, ...) {
  cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  NextMethod()
}
