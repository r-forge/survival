#SCCS $Id: untangle.specials.s,v 4.1 1992-07-10 16:36:24 therneau Exp $
untangle.specials <- function(tt, special, order=1) {
    #
    # There was a change in the specials, so action depends on your release
    #   of S
    #
    spc <- attr(tt, 'specials')[[special]]
    if (length(spc)==0)
	return(vars=character(0), terms=numeric(0))

    facs <- attr(tt, 'factor')
    fname <- dimnames(facs)

    if ((attr(terms(y~ zed(x), specials='zed'), 'specials'))$zed ==1) {
	# old style
	if (any(order>1))
	   warning("Can't select specials based on the order of the terms")
	list(vars=(fname[[2]])[spc],  terms=spc)
	}
    else {
	ff <- apply(facs[spc,,drop=F], 2, sum)
	list(vars= (fname[[1]])[spc],
	     terms= seq(ff)[ff & match(attr(tt, 'order'), order, nomatch=0)])
	}
    }
