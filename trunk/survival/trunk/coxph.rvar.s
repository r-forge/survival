#SCCS $Id: coxph.rvar.s,v 4.3 1994-01-18 10:30:48 therneau Exp $
coxph.rvar <- function(fit, collapse) {
    rcall <- match.call()
    if (class(fit) != 'coxph')
	stop ("First argument must be a fitted Cox model")

    if (missing(collapse)) temp <- residuals.coxph(fit, type='dfbeta')
    else temp <- residuals.coxph(fit, type='dfbeta', collapse=collapse)
    if (any(is.na(temp)))
       if (ncol(temp)==1) temp<- temp[!is.na(temp),,drop=F]
       else               temp <- temp[!is.na(temp %*% rep(1,ncol(temp))),]
    fit$robust.var <- t(temp) %*% temp
    fit$rcall <- rcall
    fit
    }
