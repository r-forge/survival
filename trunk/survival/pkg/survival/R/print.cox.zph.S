# $Id$
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3))
    invisible(print(x$table, digits=digits))
