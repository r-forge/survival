#SCCS $Date: 1992-04-14 18:07:40 $ $Id: residuals.coxph.null.s,v 4.2 1992-04-14 18:07:40 grill Exp $
residuals.coxph.null <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    ...)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') NextMethod()
    else stop(paste("\'", type, "\' residuals are not defined for a null model",
			sep=""))
    }
