#SCCS $Date: 1992-03-04 16:48:25 $ $Id: residuals.coxph.null.s,v 4.1 1992-03-04 16:48:25 therneau Exp $
residuals.coxreg.null <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    ...)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') NextMethod()
    else stop(paste("\'", type, "\' residuals are not defined for a null model",
			sep=""))
    }
