# $Id: frailty.s,v 1.3 2006-08-28 13:51:33 m015733 Exp $
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

			  
			   
			   
