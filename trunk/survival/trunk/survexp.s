#SCCS $Date: 1992-03-08 20:19:22 $ $Id: survexp.s,v 4.2 1992-03-08 20:19:22 therneau Exp $
surv.exp <- function(entry, birth, sex, times=round(182.6 * 0:8),
		      type=c("mean", "individual", "matrix"),
		      expected=surv.exp.uswhite, interp=F) {
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

    nomiss <- !(is.na(entry) | is.na(birth) | is.na(sex))
    if (type=='indivdual') {
	nomiss <- nomiss & !is.na(times)
	times <- times[nomiss]
	}
    else if (any(is.na(times))) stop("Missing values not allowed in 'times'")

    # toss out missings
    entry <- entry[nomiss]
    birth <-birth[nomiss]   #I need to watch out that the 'class' isn't lost
    sex <- sex[nomiss]
    nused <- sum(nomiss)
    if (nused==0) stop ("No data remains after deleting missing values")
    if (any(nomiss)) attr(nused, 'omit') <- seq(nn)[!nomiss]

    if (!inherits(birth, "date")) {
	# Assume that they gave an age
	if (any(birth <=0)) stop ("Age must be >0 ")
	if (any(birth >130))stop ("Impossibly large value for age at entry")
	birth <- entry - round(birth*365.25)
	}

    if (any(times<0)) stop ("'Times' must be >=0")
    if (type != 'individual') times <- as.integer(sort(unique(times)))
    ntime <- length(times)

    dn <- dimnames(expected)
    dd <- dim(expected)
    if (is.null(dd) || length(dd) <2 || dd[2]<2 || length(dd) >3 )
	stop("The hazard table is not in the correct format")
    ages <- as.integer(round(365.25 *as.numeric(dn[[1]])))
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
			    as.integer(dd[1]),
			    ages,
			    as.integer(dd[3]),
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

    if (any(!nomiss)) {
	omit <- seq(along=nused)[!nused]
	attr(nuse, 'omit') <- omit
	}
    if (type != 'matrix') xx <- list(time=times, surv=temp$surv, n=nused)
    else     xx <-  list(time=times,
			 surv=matrix(temp$surv, nrow=nused, byrow=T),
			 n=used)
    xx$call <- call
    if (type == 'individual') {
	if (any(!nomiss)) {
	    xx$time <- na.expand(times, omit)
	    xx$surv <- na.expand(xx$surv, omit)
	    }
	attr(xx, "class") <- "surv.exp"
	}
    else attr(xx, "class") <- c("surv.exp", "surv.fit")
    xx
    }

