#SCCS $Id: survexp.s,v 4.9 1992-04-20 17:14:33 therneau Exp $
survexp <- function(entry, birth, sex,
		      times=round(182.6 * 0:8),
		      data=sys.parent(), subset, na.action,
		      type=c("mean", "individual", "matrix"),
		      expected=survexp.uswhite, interp=F) {
    call <- match.call()

    type <- match.arg(type)

    #
    # Mimic the processing of a formula.  It is tempting to just build a
    #  formula and then use it, but then constants in the arg list don't
    #  come through.  And if you paste together I(entry) + ..., then the
    #  date attribute gets lost.
    # So evaluate, and then do the formula.  I could do the na.action and
    #  etc "by hand", but then miss out on global changes to S
    #
    if (missing(entry)) stop("The entry argument is required")
    if (missing(birth)) stop("The birth argument is required")
    if (missing(sex))   stop("The sex argument is required")
    entry <- eval(call[['entry']], data)
    birth <- eval(call[['birth']], data)
    sex   <- eval(call[['sex']], data)
    nused <- length(sex)
    if (length(birth) != nused || length(entry) != nused)
	stop("First 3 arguments must be the same length")

    m <- call
    if (type=='individual') {
	times <- eval(call[["times"]], data)
	if (length(times) != nused) stop("Wrong length for \"times\"")
	m$formula <- formula("~entry + birth + sex + times")
	}
    else {
	if (any(is.na(times))) stop("Missing values not allowed in 'times'")
	m$formula <- formula("~ entry + birth + sex")
	}
    m[[1]] <- as.name("model.frame")
    m <- m[match(c("", "formula", "subset", "na.action"),
		     names(m), 0)]
    m <- eval(m)

    entry <- m[[1]]
    birth <- m[[2]]
    sex   <- m[[3]]
    if (type=='individual') times <- m[[4]]
    else    times <- as.integer(sort(unique(times)))
    if (!inherits(entry, "date")) stop ("'Entry' must be a date")
    if (!inherits(birth, "date")) {
	# Assume that they gave an age
	if (any(birth <=0)) stop ("Age must be >0 ")
	if (any(birth >130))stop ("Impossibly large value for age at entry")
	birth <- entry - round(birth*365.2425)
	}

    if (any(entry<birth)) stop("Subject entered before birth")
    if (any(times<0)) stop ("'Times' must be >=0")
    ntime <- length(times)

    dn <- dimnames(expected)
    dd <- dim(expected)
    if (is.null(dd) || length(dd) <2 || dd[2]<2 || length(dd) >3 )
	stop("The hazard table is not in the correct format")
    ages <- as.integer(round(365.2425 *as.numeric(dn[[1]])))
    if (ages[1] !=0  ||  length(ages)<2 || any(diff(ages)<=0))
	stop("The hazard table is not in the correct format")

    if (!all(sex==floor(sex) & sex>0 & sex<=dd[2]))
	 stop("Invalid value for \'sex\'")

    if (length(dd) ==2 || dd[3]==1) {
	years <- 1
	dd[3] <- 1
	}
    else {
	years <- as.integer(dn[[3]])
	if (any(diff(years)<=0)) stop ("The hazard table must have increasing `years'")
	if (!interp || all(diff(years)==1)) years <- mdy.date(1,1, years)
	else { # linearly interpolate the rows
	    nyear <- 1+max(years) - min(years)
	    temp <- matrix(double(dd[3] * nyear), ncol=nyear)
	    diffs <- 1/diff(years)
	    new <- (min(years)):(max(years))
	    for (i in 1:(dd[3]-1)) {
		xx <- 1 - (new - years[i])*diffs[i]
		xx <- ifelse(xx<0 | xx>1, 0, xx)
		temp[i,] <- temp[i,] + xx
		temp[i+1,] <- temp[i+1,] + ifelse(xx==0, 0, 1-xx)
		}
	    temp[dd[3], nyear] <- 1
	    expected <- matrix(expected,ncol=dd[3]) %*% temp
	    dd[3] <- nyear
	    dim(expected) <- dd
	    years <- mdy.date(1,1,new)
	    }
	}
    if (type != 'matrix') nsurv <- 1
    else                  nsurv <- nused
    temp <-  .C("survexp", as.integer(ntime),
			    as.integer(times),
			    as.integer(dd),
			    ages,
			    as.integer(years),
			    as.single(expected),
			    as.integer(nused),
			    as.integer(entry),
			    as.integer(birth),
			    as.integer(sex),
			    double(ntime),
			    as.integer(nsurv),
			    as.integer(type=='individual'),
			    surv= double(ntime*nsurv))

    if (type != 'matrix') xx <- list(time=times, surv=temp$surv, n=nused)
    else     xx <-  list(time=times,
			 surv=matrix(temp$surv, ncol=nused, byrow=F),
			 n=nused)
    xx$call <- call
    if(!is.null(tj <- attr(m, 'na.action')))
	xx$na.action <- tj
    if (type == 'individual') {
	attr(xx, "class") <- "survexp"
	if (!is.null(tj)) {
	    xx$time <- naresid(tj, times)
	    xx$surv <- naresid(ty, xx$surv)
	    }
	}
    else attr(xx, "class") <- c("survexp", "survfit")
    xx
    }

