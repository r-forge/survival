#SCCS $Id: Surv.s,v 5.6 2001-06-12 08:53:57 therneau Exp $
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {

    if (missing(time)) stop ("Must have a time argument")
    if (!is.numeric(time)) stop ("Time variable is not numeric")
    nn <- length(time)
    if (length(origin)!=1 && length(origin)!=nn) 
	    stop("Wrong length for origin")

    ng <- 1 + (!missing(time2)) + (!missing(event))  #of first 3 present
    # The logic below uses "ng" throughout; why not use "missing(time2)"
    # and missing(status) instead?  Because we want to assume that 
    # "Surv(a,b)" has the variable b matched to status rather than time2.
    #

    if (missing(type)) {
	if (ng==1 || ng==2) type <- 'right'
	else if (ng==3)     type <- 'counting'
	}
    else {
	type <- match.arg(type)
	if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
	if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
		stop("Wrong number of args for this type of survival data")
	}

    if (ng==1) {
	ss <- cbind(time=time-origin, status=1)
	}
    else if (type=='right' || type=='left') {
	if (length(time2) != nn) stop ("Time and status are different lengths")
	if (is.logical(time2)) status <- 1*time2
	else  if (is.numeric(time2)) {
	    who2 <- !is.na(time2)
	    if (max(time2[who2]) ==2) status <- time2 -1
	    else status <- time2
	    if (any(status[who2] !=0  & status[who2]!=1))
				stop ("Invalid status value")
	    }
	else stop("Invalid status value")
	ss <- cbind(time=time-origin, status=status)
	}
    else  if (type=='counting') {
	if (length(time2) !=nn) stop ("Start and stop are different lengths")
	if (length(event)!=nn) stop ("Start and event are different lengths")
	if (!is.numeric(time2)) stop("Stop time is not numeric")
	who3 <- !(is.na(time) | is.na(time2))
	if (any (time[who3]>= time2[who3]))stop("Stop time must be > start time")
	if (is.logical(event)) status <- 1*event
	    else  if (is.numeric(event)) {
		who2 <- !is.na(event)
		if (max(event[who2])==2) status <- event - 1
		else status <- event
		if (any(status[who2] !=0  & status[who2]!=1))
				stop("Invalid status value")
		}
	    else stop("Invalid status value")
	ss <- cbind(start=time-origin, stop=time2-origin, status)
	}

    else {  #interval censored data
	if (type=='interval2') {
	    if (!is.numeric(time2)) stop("Time2 must be numeric")
	    event <- ifelse(is.na(time), 2,
		     ifelse(is.na(time2),0,
		     ifelse(time==time2, 1,3)))
	    if (any(time[event==3] > time2[event==3]))
		stop("Invalid interval: start > stop")
	    time <- ifelse(event!=2, time, time2)
	    type <- 'interval'
	    }
	else {
	    temp <- event[!is.na(event)]
	    if (!is.numeric(temp)) stop("Status indicator must be numeric")
	    if (length(temp)>0 && any(temp!= floor(temp) | temp<0 | temp>3))
		stop("Status indicator must be 0, 1, 2 or 3")
	    }
	status <- event
	ss <- cbind(time1=time-origin, 
		    time2=ifelse(!is.na(event) & event==3, time2-origin, 0),
		    status=status)
	}

    attr(ss, "type")  <- type
    oldClass(ss) <- 'Surv'
    ss
    }

print.Surv <- function(x, quote=F, ...) {
    print(as.character.Surv(x), quote=quote, ...)
    invisible(x)
    }

as.character.Surv <- function(xx) {
    oldClass(xx) <- NULL
    type <- attr(xx, 'type')
    if (type=='right') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste(format(xx[,1]), temp, sep='')
	}
    else if (type=='counting') {
	temp <- xx[,3]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	paste('(', format(xx[,1]), ',', format(xx[,2]), temp,
			 ']', sep='')
	}
    else if (type=='left') {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "<"," "))
	paste(temp, format(xx[,1]), sep='')
	}
    else {   #interval type
	stat <- xx[,3]
	temp <- c("+", "", "-", "]")[stat+1]
	temp2 <- ifelse(stat==3,
			 paste("[", format(xx[,1]), ", ",format(xx[,2]), sep=''),
			 format(xx[,1]))
	ifelse(is.na(stat), "NA", paste(temp2, temp, sep=''))
	}
    }

"[.Surv" <- function(x, ..., drop=F) {
    # If only 1 subscript is given, the result will still be a Surv object
    #  If the second is given extract the relevant columns as a matrix
    if (missing(..2)) {
	temp <- oldClass(x)
	type <- attr(x, "type")
	oldClass(x) <- NULL
	x <- x[..., drop=F]
	oldClass(x) <- temp
	attr(x, "type") <- type
	x
	}
    else {
	oldClass(x) <- NULL
	NextMethod("[")
	}
    }

is.na.Surv <- function(x) {
    as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")
is.Surv <- function(x) inherits(x, 'Surv')
