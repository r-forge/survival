# SCCS $Id: coxph.control.s,v 5.2 1999-06-18 15:50:06 therneau Exp $
#
# Gather all of the control parameters for coxph into one spot
#
coxph.control <- function(eps=1e-4, toler.chol =eps/10000, iter.max=10,
			  toler.inf= sqrt(eps), outer.max=10 ) {
    if (iter.max <0) stop("Invalid value for iterations")
    if (eps <=0) stop ("Invalid convergence criteria")
    if (eps <= tol.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (toler.inf <=0) stop ("The inf.warn setting must be >0")
    list(eps=eps, toler.chol=toler.chol, iter.max=iter.max, 
	 toler.inf=toler.inf, outer.max=outer.max)
    }
