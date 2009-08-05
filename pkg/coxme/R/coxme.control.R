# $Id: coxme.control.s,v 1.6 2004/03/23 19:59:51 therneau Exp $
#
# Gather all of the control parameters for coxme into one spot
#
coxme.control <- function(eps=1e-8, 
                          toler.chol = .Machine$double.eps ^ .75, 
			  inner.iter=5,
			  sparse.calc=NULL,
                          optpar=list(method='BFGS'),
                          toler.opt=1e-3 )  {
    if (iter.max <0) stop("Invalid value for iterations")
    if (inner.iter<1) stop("Invalid value for inner iterations")
    if (eps <=0) stop ("Invalid convergence criteria")
    if (eps <= toler.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (toler.opt <=eps)
        warning("For numerical accuracy, eps should be < toler.opt")
    if (!is.null(sparse.calc)) {
        if (sparse.calc !=0 && sparse.calc !=1)
            stop("Invalid value for sparse.calc option")
        }
    list(eps=eps, toler.chol=toler.chol, iter.max=iter.max,
	 inner.iter=inner.iter, sparse.calc=sparse.calc,
         optpar=c(optpar, control=list(reltol=toler.opt))
	 )
    }
