#SCCS $Date: 1998-09-25 22:43:59 $ $Id: Surv.s,v 5.2 1998-09-25 22:43:59 therneau Exp $
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {
    nn <- length(time)
    ng <- nargs()
    if (missing(type)) {
	if (ng==1 || ng==2) type <- 'right'
	else if (ng==3)     type <- 'counting'
	else stop("Invalid number of arguments")
	}
    else {
	type <- match.arg(type)
	ng <- ng-1
	if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
	if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
		stop("Wrong number of args for this type of survival data")
	}
    who <- !is.na(time)

    if (ng==1) {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	ss <- cbind(time, 1)
	dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else if (type=='right' || type=='left') {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
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
	 ss <- cbind(time, status)
	 dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else  if (type=='counting') {
	if (length(time2) !=nn) stop ("Start and stop are different lengths")
	if (length(event)!=nn) stop ("Start and event are different lengths")
	if (!is.numeric(time))stop("Start time is not numeric")
	if (!is.numeric(time2)) stop("Stop time is not numeric")
	who3 <- who & !is.na(time2)
	if (any (time[who3]>= time2[who3]))stop("Stop time must be > start time")
	if (is.logical(event)) status <- 1*event
	    else  if (is.numeric(event)) {
		who2 <- !is.na(event)
		status <- event - min(event[who2])
		if (all(status==0)) status <- status +1
		if (any(status[who2] !=0  & status[who2]!=1))
				stop("Invalid status value")
		}
	    else stop("Invalid status value")
	ss <- cbind(time-origin, time2-origin, status)
	}

    else {  #interval censored data
	if (type=='interval2') {
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
	ss <- cbind(time, ifelse(!is.na(event) & event==3, time2, 1),
			    status)
	}

    attr(ss, "type")  <- type
    oldClass(ss) <- 'Surv'
    ss
    }

print.Surv <- function(xx, quote=F, ...)
    invisible(print(as.character.Surv(xx), quote=quote, ...))

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
