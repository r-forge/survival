# SCCS $Id: print.cox.zph.s,v 4.3 1993-10-01 15:55:49 therneau Exp $
print.cox.zph <- function(x, digits = .Options$digits - 3)
    invisible(print(x$table, digits=digits))
