#SCCS $Date: 1992-03-24 22:08:43 $ $Id: predict.coxph.s,v 4.2 1992-03-24 22:08:43 therneau Exp $
predict.coxreg <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels(object), miss.expand=T, collapse, ...)

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
