#SCCS $Id: survreg.control.s,v 4.2 1998-10-27 17:55:12 therneau Exp $
survreg.control <- function(maxiter=30, rel.tolerance=1e-5, failure=1,
			    toler.chol=1e-10)
    list(maxiter = maxiter, rel.tolerance = rel.tolerance, 
	 failure =failure, toler.chol= toler.chol)
