# SCCS $Id: frailty.gaussian.s,v 1.3 1999-01-14 09:40:36 therneau Exp $
# 
# Defining function for gaussian frailty fits
#
frailty.gaussian <- function(x, sparse=(nclass >5), theta, df, 
		   method=c("reml", "aic", "df", "fixed"), ...) {

    nclass <- length(unique(x))
    # Check for consistency of the arguments
    if (missing(method)) {
	if (!missing(theta)) {
	    method <- 'fixed'
	    if (!missing(df)) 
		    stop("Cannot give both a df and theta argument")
	    }
	else if (!missing(df)) {
	    if (df==0) method <- "aic"
	    else       method <- 'df'
	    }
	}
    method <- match.arg(method)
    if (method=='df' && missing(df)) stop("Method = df but no df argument")
    if (method=='fixed' && missing(theta))
	    stop("Method= fixed but no theta argument")
    if (method !='fixed' && !missing(theta)) 
	    stop("Method is not 'fixed', but have a theta argument")

    if (sparse){
	x <-as.numeric(as.factor(x))
	oldClass(x) <- "coxph.penalty"
	}
    else{
	x <- as.factor(x)
	oldClass(x) <- "coxph.penalty"
	attr(x,'contrasts') <- function(n,...) contr.treatment(n,F)
        }
    if (!missing(theta) & !missing(df)) 
	    stop("Cannot give both a df and theta argument")

    pfun<- function(coef, theta, ndead){
	if (theta==0) list(recenter=0, penalty=0, flag=T)
	  else {
	      recenter <- mean(coef)
	      coef <- coef - recenter
	      list(recenter = recenter,
		   first=   coef/theta,
		   second=  rep(1, length(coef))/theta,
#		   penalty= -sum(log(dnorm(coef,0, sqrt(theta))),
                   penalty= 0.5* sum(coef^2/theta + log(2*pi*theta)),
		   flag=F)
	           }
	   }

     printfun <- function(coef, var, var2, df, history) {
	if (!is.null(history$history)) 
	     theta <- history$history[nrow(history$history),1]
	else theta <- history$theta
		
	if (is.matrix(var)) test <- coxph.wtest(var, coef)$test
	else 		    test <- sum(coef^2/var)
	df2 <- max(df, .5)      # Stop silly p-values
	list(coef=c(NA, NA, NA, test, df, 1-pchisq(test, df2)),
		 history=paste("Variance of random effect=", format(theta)))
	}

   if (method=='reml') {
	temp <- list(pfun=pfun, 
		     printfun=printfun,
		     diag =T,
		     sparse= sparse,
		     cargs = c('coef', 'trH', 'loglik'),
		     cfun = frailty.controlgauss,
		     cparm= list( ...))
	}
    else if (method=='fixed') {
	temp <- list(pfun=pfun,
		     printfun = printfun,
		     diag =T,
		     sparse= sparse,
		     cfun = function(parms, iter, old){
		          list(theta=parms$theta, done=T)},
		     cparm= list(theta=theta, ...))
        }
    else if (method=='aic') {
	temp <- list(pfun=pfun,
		     printfun=printfun,
		     diag =T,
		     sparse= sparse,
		     cargs = c("neff", "df", "plik"),	
		     cparm=list(...),
		     cfun = frailty.controlaic)
	}
    else {  #df method
	temp <- list(pfun=pfun,
		     printfun =printfun,
		     diag =T,
		     sparse= sparse,
		     cargs=('df'),
		     cparm=list(df=df, thetas=0, dfs=0,
		                guess=3*df/length(unclass(x)), ...),
                     cfun = frailty.controldf)
	}

    # If not sparse, give shorter names to the coefficients, so that any
    #   printout of them is readable.
    if (!sparse) {
	vname <- paste("gamma", levels(x), sep=':')
	temp <- c(temp, list(varname=vname))
	}
    attributes(x) <- c(attributes(x), temp)
    x
    }

			  
			   
			   
