#SCCS $Id: coxph.rvar.s,v 4.2 1994-01-17 09:36:32 therneau Exp $
coxph.rvar <- function(fit, collapse) {
    rcall <- match.call()
    if (class(fit) != 'coxph')
	stop ("First argument must be a fitted Cox model")

    if (missing(collapse)) temp <- residuals.coxph(fit, type='dfbeta')
    else temp <- residuals.coxph(fit, type='dbeta', collapse=collapse)

    fit$robust.var <- t(temp) %*% temp
    fit$rcall <- rcall
    fit
    }
