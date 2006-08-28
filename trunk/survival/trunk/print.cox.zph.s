# $Id: print.cox.zph.s,v 4.6 2006-08-28 14:24:11 m015733 Exp $
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3))
    invisible(print(x$table, digits=digits))
