#SCCS $Date: 1992-04-14 18:07:15 $ $Id: predict.coxph.s,v 4.4 1992-04-14 18:07:15 grill Exp $
predict.coxph <-
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
