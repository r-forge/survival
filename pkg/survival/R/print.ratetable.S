## $Id$

print.ratetable <- function(x, ...) {
  cat ("Rate table with dimension(s):", attr(x, 'dimid'), "\n")
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  NextMethod()
}
