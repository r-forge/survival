#SCCS $Id: untangle.specials.s,v 4.3 1998-08-30 16:04:24 therneau Exp $
untangle.specials <- function(tt, special, order=1) {

    spc <- attr(tt, 'specials')[[special]]
    if (length(spc)==0)
	return(list(vars=character(0), terms=numeric(0)))

    facs <- attr(tt, 'factor')
    fname <- dimnames(facs)
    ff <- apply(facs[spc,,drop=F], 2, sum)
    list(vars= (fname[[1]])[spc],
	     terms= seq(ff)[ff & match(attr(tt, 'order'), order, nomatch=0)])
    }
