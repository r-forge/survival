# SCCS $Id: frailty.s,v 1.1 1998-10-28 08:52:00 therneau Exp $
# 
# Parent function for frailty, calls the actuall working functions
#
frailty <- function(x, distribution = 'gamma', ...) {
    dlist <- c("gamma", "gaussian", "cauchy")
    i <- pmatch(distribution, dlist)
    if (!is.na(i)) distribution <- dlist[i]

    temp <- paste("frailty", distribution, sep='.')
    if (!exists(temp))
	    stop(paste("Function '", temp, "' not found", sep=""))
    (get(temp))(x, ...)
    }

			  
			   
			   
