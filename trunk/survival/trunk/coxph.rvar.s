#SCCS $Id: coxph.rvar.s,v 4.1 1993-07-04 16:08:15 therneau Exp $
coxph.rvar <- function(fit, collapse) {
    rcall <- match.call()
    if (class(fit) != 'coxph')
	stop ("First argument must be a fitted Cox model")

    if (missing(collapse)) temp <- residuals.coxph(fit, type='dbeta')
    else temp <- residuals.coxph(fit, type='dbeta', collapse=collapse)

    fit$robust.var <- t(temp) %*% temp
    fit$rcall <- rcall
    fit
    }
