# SCCS $Id: print.cox.zph.s,v 4.5 1996-09-27 11:01:55 boos Exp $
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3))
    invisible(print(x$table, digits=digits))
