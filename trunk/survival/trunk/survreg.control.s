#SCCS $Id: survreg.control.s,v 4.1 1992-11-19 11:09:48 therneau Exp $
survreg.control <- function(maxiter=30, rel.tolerance=1e-5, failure=1)
    list(maxiter = maxiter, rel.tolerance = rel.tolerance, failure =failure)
