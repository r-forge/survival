# SCCS $Id: frailty.s,v 1.2 1999-01-14 09:38:59 therneau Exp $
# 
# Parent function for frailty, calls the actuall working functions
#
frailty <- function(x, distribution = 'gamma', ...) {
    dlist <- c("gamma", "gaussian", "t")
    i <- pmatch(distribution, dlist)
    if (!is.na(i)) distribution <- dlist[i]

    temp <- paste("frailty", distribution, sep='.')
    if (!exists(temp))
	    stop(paste("Function '", temp, "' not found", sep=""))
    (get(temp))(x, ...)
    }

			  
			   
			   
