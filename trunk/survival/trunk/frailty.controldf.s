# SCCS $Id: frailty.controldf.s,v 1.1 1998-10-28 08:51:56 therneau Exp $
# A function to calibrate the df
#    very empirical  
# Find the closest 3 points that span the target value
#   We know the function is monotone, so fit the function
#     dy = a * (dx)^p   to the 3 points, where dx and dy are the distance
#   from the leftmost of the three points.
# On input, parms$df = target degrees of freedom
#           parms$dfs, parms$thetas = known values (usually 0,0)
#           parms$guess = first guess
#
frailty.controldf <- function(parms, iter, old, df) {
    if (iter==0) {  
	theta <- parms$guess
	return(list(theta=theta, done=F, 
		    history=cbind(thetas=parms$thetas, dfs=parms$dfs)))
	}

    eps <- parms$eps
    if (length(eps)==0) eps <- .1

    thetas <- c(old$history[,1], old$theta)
    dfs    <- c(old$history[,2], df)
    nx <- length(thetas)
    if (nx==2) {
	#linear guess based on first two 
	# but try extra hard to bracket the root
	theta <- thetas[1] + (thetas[2]-thetas[1])*(parms$df - dfs[1])/
						    (dfs[2] - dfs[1])
	if (parms$df > df) theta <- theta * 1.5 
	return(list(theta=theta, done=F,
		    history=cbind(thetas=thetas, dfs=dfs)))
	}
    else{
	# Now, thetas= our guesses at theta
	#  dfs = the degrees of freedom for each guess
	done <- (iter>1 &&
		 (abs(dfs[nx]-parms$df) < eps))

	# actually look for a new minimum
	x <- thetas
	y <- dfs
	target <- parms$df
	if ((x[1]-x[2])*(y[1]-y[2]) >0)  y <- sort(y)  #monotone up
	else  { #monotone down
	    y <- sort(-y)
	    target <- -target
	    }
	x <- sort(x)

	if (all(y>target)) b1 <- 1     #points 1:3 are the closest then
	else if (all(y<target)) b1 <- nx-2
	else {
	    b1 <- max((1:nx)[y <= target]) #this point below target, next above
	    # use either b1,b1+1,b1+2 or  b1-1, b1, b1+1, whichever is better
	    #  better = midpoint of interval close to the target

	    if ((b1+1)==nx ||
		(b1>1 &&  ((target -y[b1]) < (y[b1+1] -target))))
		    b1 <- b1-1
	    }

	#now have the best 3 points
	# fit them with a power curve anchored at the leftmost one
	b2 <- b1 + 1:2	
	xx <- log(x[b2] - x[b1])
	yy <- log(y[b2] - y[b1])
	power <- diff(yy)/diff(xx)
	a <- yy[1] - power*xx[1]
	newx <- (log(target -y[b1]) - a)/power
	if (length(parms$trace) && parms$trace){
	    print(cbind(thetas=thetas, dfs=dfs))
	    cat("  new theta=" , format(newx), "\n\n")
	    }
	list(theta=x[b1] + exp(newx), done=done, 
	     history=cbind(thetas=thetas, dfs=dfs))
	}
    }

