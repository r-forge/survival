#SCCS $Id: survexp.s,v 4.10 1992-05-28 21:46:33 therneau Exp $
survexp <- function(entry, birth, sex,
		      times=round(182.6 * 0:8),
		      data=sys.parent(), subset, na.action,
		      type=c("mean", "individual", "matrix"),
		      expected=survexp.uswhite, interp=F) {
    call <- match.call()

    type <- match.arg(type)

    m <- call
    if (type=='individual') {
	if (missing(times))
	    stop ("The 'times' argument is required for individual survival")
	temp <- c("entry", "birth", "sex", "times")
	}
    else {
	temp <- c("entry", "birth", "sex")
	if (any(is.na(times))) stop("Missing values not allowed in 'times'")
	}
    temp2 <- NULL
    for (i in temp) temp2 <- c(temp2, deparse(call[[i]]))
    m$formula <- formula(paste("~", paste(temp2, collapse="+")))
    m[[1]] <- as.name("model.frame")
    m <- m[match(c("", "formula", "subset", "na.action", "data"),
		     names(m), 0)]
    m <- eval(m)

    entry <- m[[1]]
    birth <- m[[2]]
    sex   <- m[[3]]
    if (type=='individual') times <- m[[4]]
    else    times <- as.integer(sort(unique(times)))
    n <- length(entry)
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
    else                  nsurv <- n
    temp <-  .C("survexp", as.integer(ntime),
			    as.integer(times),
			    as.integer(dd),
			    ages,
			    as.integer(years),
			    as.single(expected),
			    as.integer(n),
			    as.integer(entry),
			    as.integer(birth),
			    as.integer(sex),
			    double(ntime),
			    as.integer(nsurv),
			    as.integer(type=='individual'),
			    surv= double(ntime*nsurv))

    if (type != 'matrix') xx <- list(time=times, surv=temp$surv, n=n)
    else     xx <-  list(time=times,
			 surv=matrix(temp$surv, ncol=n, byrow=F),
			 n=n)
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

