.BG
.FN survsum
.TL
Kaplan-Meier survival curve and percentages at selected times
.DN
Calculates Kaplan-Meier survival percentages, standard error, and number 
at risk at specified times for defined groups.  Test for difference between
survival curves.
.CS
survsum(formula, data=sys.parent(), sptms=NULL, xlim, 
tlines=T, log=F, xscale=1,yscale=100, mark.time=F, mark=3, 
cex=1, xlab="Time", ylab="Survival (%)", lgd="bl",
ttl="",...)
.RA
.AG formula
a formula expression as for other survival models, of the  form
Surv(time,  status) ~ predictors.  Same formula expression as used 
in 'survfit'. See help file for survfit.  Maximum of 6 groups.
.OA
.AG data
a data.frame in which to interpret the variables named in the formula.
.AG sptms
a specified vector of positive times at which to compute the survival 
percentages, standard errors, and numbers at risk.  A maximum of four 
times can be specified.
.AG xlim
a vector of the form: c(x1,x2).  The approximate minimum and maximum values 
to  be  put  on  x-axis.  Default sets x1=0 and x2=maximum time value.
.AG tlines
a logical value indicating whether vertical lines and labels should be drawn
on the plot at the specified times.
.AG log
logical value: should the y axis be on a log scale?
.AG xscale
a scalar to be used to divide the x axis.  A value of 365, for instance, would
be used to convert from days to years.
.AG yscale
a scalar to be used to multiply the y axis.  The default value of 100 is used
to get a percent scale.  'yscale=1' would set the y axis from 0 to 1.
.AG mark.time
controls the labeling of the curves.  If set to True then curves are 
marked at each censoring time.  If mark.time is a numeric vector, then curves
are marked at these specified time points.
.AG mark
vector of mark parameters, which will be used to  label the  curves.  The
'lines' help file contains examples of the possible marks.  The 
vector is resued cyclically if it  is shorter than the number of curves.
.AG cex
parameter available to change the size of "mark".   Not a vector; all marks
have the same size.
.AG xlab
character string label for the x axis.
.AG ylab
character string label for the y axis.
.AG lgd
legend placement.  "tr"=top right corner of the plot, "under"=under the plot,
"n" omits the legend.  The default is: "bl"=bottom left corner of the plot.
.AG ttl
title to be printed in upper left corner of page (only visible on printed
copy).
.AG ...
In addition, the high-level graphics arguments described under par and the
arguments to title may be supplied to this function.
.RT
A list of class `"htest"', containing the following components:
.PP
.RC no.pts
the total number of data points in each group.
.RC no.events
the total number of events in each group.
.RC chisq
the chisquare statistic for the test of a difference between survival curves.
.RC p
the p-value for the above test.
.RC t1
a matrix containing the survival percentages, standard errors, and numbers at
risk for all groups at time t1.
.RC t2,t3,t4
see above.
.SE
A plot with multiple survival curves is drawn (one for each group). This plot
includes: overall group statistics and group statistics at each specified 
time point and a test for a difference between survival curves.
.DT
The total number of points and events is reported for each group.  For each 
specified time point, group survival percentages (followed by standard error
and number left at risk) are computed.
.PP
The test for a difference bewteen survival curves uses the chisquare 
statistic from the 'survdiff' function with rho=0.  This is the log-rank test.
.SH AUTHOR
Mark Dietrich, Mayo Clinic Section of Medical Research Statistics
summer student 1992.
.SA
Surv survdiff survexp survfit
.EX
survsum (Surv (futime,dead)~sex+dose,data=blood.Dat,sptms=c(20,30,45),
xlim=c(0,50))
##groups are all combinations of 'sex' and 'dose', specified times are 20, 30,
and 45.  
.KW survival

.WR
