#SCCS $Date: 1992-03-25 01:33:56 $ $Id: Surv.s,v 4.4 1992-03-25 01:33:56 therneau Exp $
# Package up surivival type data as a structure
#  Eventually allow lots of censored data types
#
Surv <- function(time, time2, event,
		       type=c('right', 'left', 'interval', 'counting')) {
    nn <- length(time)
    ng <- nargs()
    if (missing(type)) {
	if (ng<3) type <- 'right'
	else      type <- 'counting'
	}
    else {
	type <- match.arg(type)
	ng <- ng-1
	}
    who <- !is.na(time)

    if (ng==1) {
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	else if (any(time[who]<0))  stop ("Time variable must be >= 0")
	ss <- cbind(time, 1)
	dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else if (ng==2) {  #assume second arg is event
	if (!is.numeric(time)) stop ("Time variable is not numeric")
	else if (any(time[who]<0))  stop ("Time variable must be >= 0")
	if (length(time2) != nn) stop ("Time and status are different lengths")
	if (is.logical(time2)) status <- 1*time2
	    else  if (is.numeric(time2)) {
		who2 <- !is.na(time2)
		status <- time2 - min(time2[who2])
		if (any(status[who2] !=0  & status[who2]!=1))
				stop ("Invalid status value")
		}
	    else stop("Invalid status value")
	 ss <- cbind(time, status)
	 dimnames(ss) <- list(NULL, c("time", "status"))
	}
    else  {
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
		if (any(status[who2] !=0  & status[who2]!=1))
				event ("Invalid status value")
		}
	    else event("Invalid status value")
	ss <- cbind(time, time2,status)
	}
    attr(ss, "class") <- c("Surv")
    attr(ss, "type")  <- type
    ss
    }

print.Surv <- function(xx, quote=F, ...) {
    class(xx) <- NULL
    if (ncol(xx)==2) {
	temp <- xx[,2]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	print(paste(xx[,1], temp, sep=''), quote=quote)
	}
    else {
	temp <- xx[,3]
	temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
	print(paste('(', xx[,1], ',', xx[,2], temp, ']', sep=''), quote=quote)
	}
    }

"[.Surv" <- function(x,i,j, drop=F) {
    temp <- class(x)
    type <- attr(x, "type")
    class(x) <- NULL
    attr(x, 'type') <- NULL
    if (missing(j)) {
	x <- x[i,,drop=drop]
	class(x) <- temp
	attr(x, "type") <- type
	x
	}
    else NextMethod("[")
    }

is.na.Surv <- function(x) {
    class(x) <- NULL
    as.vector( (1* is.na(x))%*% rep(1, ncol(x)) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")