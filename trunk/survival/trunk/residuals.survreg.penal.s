# SCCS $Id: residuals.survreg.penal.s,v 1.1 1998-11-30 08:23:06 therneau Exp $
# This routine just stops disastrous arithmetic for models with sparse
# terms.  A placeholder until the proper sparse terms actions are inserted.
residuals.survreg.penal <- function(object, ...) {
    pterms <- object$pterms
    if (any(pterms==2))
	    stop("Residualss not available for sparse models")
    NextMethod('residuals')
    }
