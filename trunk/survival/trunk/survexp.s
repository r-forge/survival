#SCCS $Id: survexp.s,v 4.6 1992-03-30 13:21:02 therneau Exp $
surv.exp <- function(entry, birth, sex, times=round(182.6 * 0:8),
		      type=c("mean", "individual", "matrix"),
		      expected=surv.exp.uswhite, interp=F,
		      na.action= na.omit) {
    call <- match.call()
    if (!inherits(entry, "date")) stop ("'Entry' must be a date")
    nn <- length(entry)
    if (length(birth) != nn) stop("Entry and birth not the same length")
    if (length(sex) != nn) stop ("Entry and sex not the same length")

    type <- match.arg(type)
    if (missing(type)){
	if (!missing(times) && length(times)==nn) type <-'individual'
	else type <- 'mean'
	}

    # Use na.action to deal with missing values
    if (missing(na.action) &&
	      !is.null(tj <- options("na.action")[[1]])) na.action <- get(tj)
    if (type=='individual') {
	temp <- list(entry, birth, sex, times)
	names(temp) <- names(call)[1:4]
	attr(temp, 'row.names') <- 1:nn
	attr(temp, 'class') <- "data.frame"
	m <- na.action(temp)
	times <- m[[4]]
	}
    else {
	if (any(is.na(times))) stop("Missing values not allowed in 'times'")
	temp <- list(entry, birth, sex)
	names(temp) <- names(call)[1:3]
	attr(temp, 'row.names') <- 1:nn
	attr(temp, 'class') <- "data.frame"
	m <- na.action(temp)
	}
    entry <- m[[1]]
    birth <- m[[2]]
    sex   <- m[[3]]
    nused <- length(sex)

    if (!inherits(birth, "date")) {
	# Assume that they gave an age
	if (any(birth <=0)) stop ("Age must be >0 ")
	if (any(birth >130))stop ("Impossibly large value for age at entry")
	birth <- entry - round(birth*365.2425)
	}

    if (any(entry<birth)) stop("Subject entered before birth")
    if (any(times<0)) stop ("'Times' must be >=0")
    if (type != 'individual') times <- as.integer(sort(unique(times)))
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
	attr(xx, "class") <- "surv.exp"
	if (!is.null(tj)) {
	    xx$time <- naresid(tj, times)
	    xx$surv <- naresid(ty, xx$surv)
	    }
	}
    else attr(xx, "class") <- c("surv.exp", "surv.fit")
    xx
    }

