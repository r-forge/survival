#SCCS $Date: 1992-03-04 16:48:14 $ $Id: predict.coxph.s,v 4.1 1992-03-04 16:48:14 therneau Exp $
predict.coxreg <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels(object), miss.expand=T, collapse, ...)

    {
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
