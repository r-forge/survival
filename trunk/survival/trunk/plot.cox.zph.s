#SCCS $Id: plot.cox.zph.s,v 4.1 1993-03-04 11:58:02 therneau Exp $
plot.cox.zph <- function(x, resid=T, se=T, df=4, nsmo=40) {
    d <- nrow(x$y)
    nvar <- ncol(x$y)
    pred.x <- seq(from=min(x$x), to=max(x$x), length=nsmo)
    temp <- c(pred.x, x$x)
    lmat <- ns(temp, df=df, intercept=T)
    pmat <- lmat[1:nsmo,]       # for prediction
    xmat <- lmat[-(1:nsmo),]
    qmat <- qr(xmat)

browser()
    if (se) {
	bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
	xtx <- bk %*% t(bk)
	temp <- xtx %*% ( (d+1)*t(xmat) - c(rep(1,d)%*%xmat)) %*% xmat %*% xtx
	seval <- ((pmat%*%temp) *pmat) %*% rep(1, df)
	}

    ylab <- paste("Beta(t) for", dimnames(x$y)[[2]])
    if (x$transform=='identity') xlab <- "Time"
    else                         xlab <- paste(x$transform, '(time)', sep='')

    for (i in 1:nvar) {
	y <- x$y[,i]
	yhat <- pmat %*% qr.coef(qmat, y)
	if (resid) yr <-range(yhat, y)
	else       yr <-range(yhat)
	if (se) {
	    temp <- 2* sqrt(x$var[i,i]*seval)
	    yup <- yhat + temp
	    ylow<- yhat - temp
	    yr <- range(yr, yup, ylow)
	    }

	plot(range(x$x), yr, type='n', xlab=xlab, ylab=ylab[i])
	if (resid) points(x$x, y)
	lines(pred.x, yhat)
	if (se) {
	    lines(pred.x, yup,lty=2)
	    lines(pred.x, ylow, lty=2)
	    }
	}
    }
