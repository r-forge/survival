# SCCS $Id: print.cox.zph.s,v 4.4 1996-09-27 10:35:03 boos Exp $
print.cox.zph <- function(x, digits = max(options()$digits - 4, 3))
    invisible(print(x$table, digits=digits))
