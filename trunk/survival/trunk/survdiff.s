#SCCS $Date: 1992-06-01 09:54:04 $ $Id: survdiff.s,v 4.7 1992-06-01 09:54:04 therneau Exp $
survdiff <- function(formula, data, subset, rho=0) {
    call <- match.call()
    m <- match.call(expand=F)
    m$rho <- NULL

    if (!inherits(formula,"formula")) {
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	if (inherits(formula,"Surv")) {
	    xx <- function(x) formula(x)
	    m$formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
	    }
	else stop("Invalid formula")
	}
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- attr(m, 'terms')
    y <- model.extract(m, response)
    ny <- ncol(y)
    n <- nrow(y)

    offset<- attr(Terms, "offset")
    if (!is.null(offset)) {
	#one sample test
	offset <- as.numeric(m[[offset]])
	if (length(Terms)>0) stop("Cannot have both an offset and groups")
	if (any(offset <0 | offset >1))
	    stop("The offset must be a survival probability")
	expected <- sum(1-offset)
	observed <- sum(y[,ny])
	if (rho!=0) {
	    num <- sum(1/rho - ((1/rho + y[,ny])*offset^rho))
	    var <- sum(1- offset^(2*rho))/(2*rho)
	    }
	else {
	    var <-  sum(-log(offset))
	    num <-  var - observed
	    }
	chi <- num*num/var
	}

    else {
	#k sample test
	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) strats <- factor(rep(1,n))
	else {
	    temp <-  rep(1, length(ll))
	    strat <- attr(Terms, 'specials')$strata
	    if (length(strat)) temp[strat] <- 0
	    lname <- ifelse(temp==1, paste(ll,'=', sep=''), "")
	    temp <- paste("'",lname, "', ","as.character(m[['", ll, "']])", sep='')
	    temp <- paste(temp, collapse=", ', ',")
	    temp <- paste("paste(", temp, ",sep='')")
	    strats <- factor(eval(parse(text=temp)))
	    }
	strat2 <- as.numeric(strats)
	ngroup <- max(strat2)
	if (ngroup <2) stop ("There is only 1 group")

	if (ny!=2) stop("Surf.diff does not apply to interval time data (yet)")
	ord <- order(y[,1])

	xx <- .C("survdiff", as.integer(n),
		       as.integer(ngroup),
		       as.double(rho),
		       as.double(y[ord,1]),
		       as.integer(y[ord,2]),
		       as.integer(strat2[ord]),
		       observed = double(ngroup),
		       expected = double(ngroup),
		       var.e    = double((ngroup-1) * (ngroup-1)),
		       double(ngroup), double(n))
	n <- table(strats)
	expected <- xx$expected
	observed <- xx$observed
	temp2 <- (observed - expected) [-1]
	chi  <- sum(solve(matrix(xx$var.e, ncol=ngroup-1), temp2) * temp2)
	}

    rval <-list(n= n, obs = observed, exp=expected,
		    chisq= chi)
    na.action <- attr(m, "na.action")
    if (length(na.action)) rval$na.action <- na.action
    attr(rval, "class") <- 'survdiff'
    rval
    }
