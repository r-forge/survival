# SCCS $Id: format.Surv.s,v 4.6 1998-12-22 08:54:21 therneau Exp $
#
format.Surv <- function(x, ...) format(as.character.Surv(x), ...)

# The function to "make something suitable for inclusion in a data frame"
#   was "as.data.frame.x" in versions <5, now it is "data.frameAux.x",
#   so here we have a version specific definition.

if (version$major >= 5) {
    data.frameAux.Surv <- function(x, ...) data.frameAux.AsIs(x, ...)
    }
else {
    as.data.frame.Surv <- as.data.frame.model.matrix
    }
