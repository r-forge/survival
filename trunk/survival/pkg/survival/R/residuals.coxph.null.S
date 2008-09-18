# $Date: 2006-08-28 15:24:59 $ $Id$
residuals.coxph.null <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    ...)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') NextMethod()
    else stop(paste("\'", type, "\' residuals are not defined for a null model",
			sep=""))
    }
