#
# Create the survreg.distributions object
#
survreg.distributions <- list(
'extreme' = list(
    name = "Extreme value",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2)/1.3)/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			if (any(status==3)) {
			    temp <- ifelse(status==3,(y[,2] - y[,1])/scale,1)
			    temp2 <- temp/(exp(temp)-1)
			    temp3 <- log(exp(-temp2) - exp(-temp2*exp(temp)))
			    best <- ifelse(status==1, -(1+log(scale)),
				    ifelse(status==3, temp3, 0))
			    }
			else best <- ifelse(status==1, -(1+log(scale)), 0)
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) =",
					       format(exp(parms)))
    ),

logistic = list(
    name  = "Logistic",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2)/2)/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			if (any(status==3)) {
			    temp <- ifelse(status==3,(y[,2] - y[,1])/scale,1)
			    temp <- (y[,2] - y[,1])/scale
			    temp2 <- exp(temp/2)
			    temp3 <- log((temp2-1)/(temp2+1))
			    best <- ifelse(status==1, -log(4*scale),
				    ifelse(status==3, temp3, 0))
			    }
			else best <- ifelse(status==1, -log(4*scale), 0)
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) est =",
					       format(exp(parms)))
    ),

gaussian = list(
    name  = "Gaussian",
    parms = list(scale=1),
    logconcave = T,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2))/2 ,0)
			matrix(temp, nrow=1, dimnames=list("Log(scale)", NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			temp <-  ifelse(status==3, (y[,2] - y[,1])/scale, 1)
			temp2 <- 2*(pnorm(temp/2) -.5)
			best <- ifelse(status==1, -log(sqrt(pi)*scale),
				ifelse(status==3, temp2, 0))
			2*(best-loglik)
			},
    print = function(parms, fixed)
		if (fixed)
		     paste("Dispersion (scale) fixed at", format(exp(parms)))
		else paste("Dispersion (scale) =",
					       format(exp(parms)))
    ),

t = list(
    name  = "Student-t",
    parms = list(scale=1, df=2),
    logconcave = F,
    init  = function(x,fixed, init) {
			if (!is.null(fixed$scale))
				temp<- c(log(fixed$scale), 1)
			else if(!is.null(init$scale))
				temp<- c(log(init$scale), 0)
			else    temp<- c(log(mean(x^2))/2 ,0)

			if (!is.null(fixed$df))
				temp<- c(temp, fixed$df, 1)
			else if(!is.null(init$df))
				temp<- c(temp, init$df, 0)
			else    {
			    k <- mean(x^4)/ (3*temp[1])
			    df<- (3*k + sqrt(8*k + k^2))/(16*k*(k-1))
			    temp<- c(temp, df, 0)
			    }

			matrix(temp, nrow=2, byrow=T,
			    dimnames=list(c("Log(scale)", "df"),  NULL))
			},
    deviance= function(y, parms, loglik) {
			status <- y[,ncol(y)]
			scale <- exp(parms[1])
			df <- parms[2]
			temp <-  ifelse(status==3, (y[,2] - y[,1])/scale, 1)
			temp2 <- 2*(pt(temp/2, df) -.5)
			temp3 <- lgamma((df+1)/2) -
				    (lgamma(df/2) + .5*log(pi*df*scale^2))
			best <- ifelse(status==1, temp3,
				ifelse(status==3, temp2, 0))
			2*(best-loglik)
			},
    print = function(parms, fixed) {
	    tt <- if (fixed[1])
		     paste("Dispersion (scale) fixed at", format(exp(parms[1])))
		else paste("Dispersion (scale) =",
					       format(exp(parms[1])))
	    if (fixed[2]) paste(tt, ", df fixed at", format(parms[2]))
	    else tt
	    }
    )
 )
