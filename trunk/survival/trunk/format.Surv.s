# SCCS $Id: format.Surv.s,v 4.5 1998-08-30 15:18:37 therneau Exp $
#
# These two functions operate with the newer data.frame code, found on statlib
#   (Eventually part of S, I assume)
#
format.Surv <- function(x, ...) format(as.character.Surv(x), ...)
as.data.frame.Surv <- as.data.frame.model.matrix
