# SCCS $Id: print.cox.zph.s,v 4.2 1993-04-06 15:53:00 therneau Exp $
print.cox.zph <- function(x, digits = .Options$digits - 3)
    print(x$table, digits=digits)
