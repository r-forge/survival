# $Id: survsum.s,v 1.1 1992-08-14 15:37:12 grill Exp $
# written by Mark Dietrich
survsum <- function	
		(formula,data=sys.parent(),sptms=NULL,xlim,tlines=T,log=F,
		xscale=1,yscale=100,mark.time=F,mark=3,cex=1,xlab="Time",
		ylab="Survival (%)",lgd="bl",ttl="",...) {
#
	ltms <- length (sptms)			##number of specified times
#
#------------------------------------------------------------------------------
#
	if (ltms >4){
        	stop("Maximum number of specified times is four.")}
	if (ltms > 0){						##warnings
		if( any (sptms < 0))
			stop ("specified times must be positive")}
#
#------------------------------------------------------------------------------
#################
#total statistics#
##################
#
	fit <- survfit (formula,data)			##survfit object
#
	n.pts <- summary (fit,times=0,print.it=F)$n.risk
					      ##number of points in each group
	strat <- fit$strata		      ##
	gp.name <- names(strat)	      ## the group names
	n <- length (n.pts)		      ##the number of groups
#
	if (n > 6) {					##too many groups
		stop("Maximum number of groups is 6")}
#
	code <- fit$n.event				##coded events by time
	stemp <- rep (1:n,strat)
		events <- 1:n		     ##number of events in each group 
        		for (i in (1:n))
                		events[i] <- sum (code [stemp ==i])
#
#------------------------------------------------------------------------------
###############
#survival plot#
###############
#
	par (fig = c(0,1,.10,.75))
	par (err=-1)			       ##supress out-of-bounds msgs
#
	frame ()
	if (missing(xlim)){		##conditional: no xlim specified
#
		if (log) {			##conditional: log=True
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",xscale=xscale,log=T,
	,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,cex=cex,las=1,...)
								         ##Plot
		} else {			##conditional: log=False
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",yaxs="i",
	xscale=xscale,yscale=yscale,xlab=xlab,ylab=ylab,mark.time=mark.time,
	mark=mark,cex=cex,ylim=c(0,(yscale+(yscale/10))),las=1,...)	##Plot
#
		}
	xlim <- c(0,max(ends$x))       	##supply xlim value needed below
#
	} else {			##conditional: xlim is specified
#
		if (log) {			##conditional: log=True
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",xlim=xlim,
	log=T,xscale=xscale,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,
	cex=cex,las=1,...)     						##Plot
#
		} else {			##conditional: log=False
#
	ends <- plot (fit,lty=c(1,2,7,8,3,5),xaxs="i",yaxs="i",
	xlim=xlim,xlab=xlab,ylab=ylab,mark.time=mark.time,mark=mark,cex=cex,
	xscale=xscale,yscale=yscale,ylim=c(0,(yscale+(yscale/10))),las=1,...)
									##Plot
			}}
#
#
	small <- xlim[1]			##minimum x value on plot 
	big <- xlim[2]				##maximum x value on plot
#	
	par (err=0)				##error msgs resumed
#
#------------------------------------------------------------------------------
#################
#specified times#
#################
#
	vec <- rep(NA,4)		##vector of NA's:used if ltms=0
	vec[1:ltms] <- sptms		##NA's replaced with specified times
	t1 <- vec[1]			##variables assigned timepoints
	t2 <- vec[2]
	t3 <- vec[3]
	t4 <- vec[4]
#
#-----------------------------------------------------------------------------
################
#vertical lines#
################
#
        if (tlines){				##conditional: tlines=True
		lin.ok <- vec [!is.na(vec)]	##times that are not NA
#
	if (ltms != 0){
	if (any (lin.ok < small) | any (lin.ok > big))
						##conditional:	times <
						##             xmin or > xmax
		stop ("a specified (sptms) time point falls outside the plot region.")
#
		abline (v=c(lin.ok),lty=4)		 ##vertical lines
                axis (3,at=c(lin.ok),labels=c(lin.ok))}} ##axis labels at times
#
#------------------------------------------------------------------------------
########
#legend#
########
#
	par (fig=c(0,1,0,1))
	a <- 1:100
	b <- 1:100
#
	plot (a,b,type="n",xlab="",ylab="",bty="n",axes=F)
#
	if (lgd=="bl") {		##conditional: legend at bottom-right
#
	if (n == 6) {
		legend (0,30,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 5) {
		legend (0,28,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 4) {
		legend (0,26,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 3) {
		legend (0,24,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
	if (n == 2) {
		legend (0,22,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
			} else {
		legend (0,20,legend=names (fit$strata),
							lty=c(1,2,7,8,3,5))
}}}}}
        } else {
#
	if (lgd =="tr") {		##conditional: legend at top-right
#
		legend (75,69,legend=names (fit$strata),
		lty=c(1,2,7,8,3,5))
#
        } else {
#
	if (lgd == "n") 			##conditional: no legend 
		{	}}}
#
	par (fig = c(0,1,0,1))
#
	if (lgd == "under"){		##conditional: legend under plot
		legend (75,4,legend=names (fit$strata),
		lty=c(1,2,7,8,3,5))}
#
#------------------------------------------------------------------------------
#####################
#test for difference#
#####################
#
	sdif <- survdiff(eval(fit$call$formula), eval (fit$call$data))
					##survdiff function
        chisq <- round(sdif$chisq,2)	##chisquare statistic
        df <- length (sdif$n) - 1	##degrees of freedom
        p <- round(1 - pchisq(sdif$chisq, df),4)	##p-value
                if (p < .0001) (p <- ".0001")		
#
	text (0,-5,paste("LOGRANK TEST (all curves equal): Chi-Square = "
	,chisq,"   ","df = ",df,"   ","p = ",p),adj=0)
#
	mtext(date(),side=1,outer=T,adj=1)		##date
#
#------------------------------------------------------------------------------
#############################
#printing on graphics device#
#############################
#
	ysfull <- c(80,72,64,56,48,40)          ##y-values
        ys <- ysfull[1:n]
#
	par (fig = c(0,1,.5,1))
#
	mtext(ttl,side=3,outer=T,adj=1)
	plot (a,b,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
	title("K-M Survival")
#
	x1 <- c(rep(-5,n+1),-5,20)
	y1 <- c(90,ys,90,105)
	labs1 <- c("Group",gp.name,"_________________________________________",
			"Totals")
	text (x1,y1,labs1,adj=0)
#
	x2 <- c(20,rep(18,n),32,rep(30,n))
	y2 <- c(90,ys,90,ys)
	labs2 <- c("No.Pts.",n.pts,"No.Events",events)
	text (x2,y2,labs2,adj=1)
#
#
##########################
#specific time statistics#
##########################
#
	gt1 <- gt2 <- gt3 <- gt4 <- NA		##return value dummy variables

	if (!is.na (t1)) {			##conditional: t1 is not NA
#
	text (38,105,"Estimated survival % (SE,n) at time (t)",adj=0)
#
#################
#define function#
#################
#
	group <- function (ti,x,current.fit,m,gpn,scale,endsx) {
#
        	        if (ti > max(endsx)){	##conditional: time > xmax
#
	        text (x,90,c(paste("t=", ti),"___________"),adj=0)
                text(x,80,"no data",adj=0)
        	        } else {		##conditional: time < xmax
#
	        mat <- summary.survfit (current.fit,times=ti,print.it=F,
					scale=scale)
#	
		gps <- mat$strata	##group names used at time

	        bigm <- matrix(rep(NA,m*3),ncol=3)
#
	        dimnames (bigm) <- list (gpn,c("percs","se","n"))
#
        	percs <- format (round (mat$surv*100,1))  ##survival percentage
        	sters <- format(round (mat$std.err*100,1))     ##standard error
        	nrisk <- format(mat$n.risk)		       ##no. at risk
#
		bigm [as.character(gps),] <- c(percs,sters,nrisk)

        	ysfull <- c(80,72,64,56,48,40)          ##y-values
        	ys <- ysfull[1:m]
#
        	text (x,90,c(paste("t=", ti),"___________"),adj=0)
        	text (rep(x,m),ys,paste(bigm[1:m,1],"(",bigm[1:m,2],",",
		bigm[1:m,3],")"),adj=0)
		list (bigm = bigm)
		}}
#
		f1 <- group (t1,35,fit,n,gp.name,xscale,ends$x)
		gt1 <- f1$bigm
		gt2 <- gt3 <- gt4 <- NA
#
	if (!is.na (t2)) {			##conditional: t2 is not NA
		f2 <- group (t2,52,fit,n,gp.name,xscale,ends$x)
		gt2 <- f2$bigm
		gt3 <- gt4 <- NA
#
	if (!is.na (t3)) {			##conditional: t3 is not NA
		f3 <- group (t3,69,fit,n,gp.name,xscale,ends$x)
		gt3 <- f3$bigm
		gt4 <- NA
#
	if (!is.na (t4)) {			##conditional: t4 is not NA
		f4 <- group (t4,86,fit,n,gp.name,xscale,ends$x)
		gt4 <- f4$bigm
}}}}
invisible(list(no.pts=n.pts,no.events=events,chisq=chisq,p=p,t1=gt1,
		t2=gt2,t3=gt3,t4=gt4 ))			##return values
}
















