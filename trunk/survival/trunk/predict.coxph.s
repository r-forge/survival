#SCCS $Date: 1992-03-30 02:40:09 $ $Id: predict.coxph.s,v 4.3 1992-03-30 02:40:09 therneau Exp $
predict.coxreg <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels(object), collapse, ...)

    {
    type <-match.arg(type)
    if (type=='terms') {
	tt <- object$terms
	sp <- attr(tt, 'specials')$strata
	if (length(sp)) {
	    tt <- tt[-sp,]
	    object$terms <- tt
	    }
	}
    NextMethod('predict')
    }
