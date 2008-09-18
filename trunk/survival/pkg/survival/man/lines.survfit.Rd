\name{lines.survfit}
\alias{lines.survfit}
\title{
  Add Lines to a Survival Plot 
}
\description{
Usually used to add the expected survival curve(s) to a Kaplan-Meier plot 
generated with \code{plot.survfit}. 
}
\usage{
lines.survfit(x, type="s", mark=3, col=1, lty=1, lwd=1, mark.time=T, 
	xscale=1,  firstx=0, firsty=1, xmax, fun, conf.int=F,  ...)
}
\arguments{
\item{x}{
a survival object, generated from the \code{survfit} or \code{survexp} functions. 
}
\item{type}{
the line type, as described in \code{lines}.  The default is a step function 
for \code{survfit} objects, and a connected line for \code{survexp} objects. 
}
\item{mark, col, lty, lwd}{
vectors giving the mark symbol, color, line type and line width for the 
added curves. 
}
\item{mark.time}{
controls the labeling of the curves.   
If \code{FALSE}, no labeling is done.   
If \code{TRUE}, then curves are marked at each censoring time.   
If \code{mark.time} is a numeric vector, then curves are marked at  
the specified time points. 
}
\item{xscale}{
a number used to divide the x values.  If time was originally in days, a 
value of 365.25 would give a plotted scale in years. 
}
\item{firstx, firsty}{
the starting point for the survival curves.  If either of these is set to 
\code{NA} or < blank > the plot will start at the first time point of the curve. 
}
\item{xmax}{
the maximum horizontal plot coordinate.     
This shortens the curve before plotting it, so unlike using the 
\code{xlim} graphical parameter, warning messages about out of bounds points are 
not generated. 
}
\item{fun}{
an arbitrary function defining a transformation of the survival curve. 
For example \code{fun=log} is an alternative way to draw a log-survival curve 
(but with the axis labeled with log(S) values). 
Four often used transformations can be specified with a character 
argument instead: "log" is the same as using the \code{log=T} option, 
"event" plots cumulative events (f(y) = 1-y), 
"cumhaz" plots the cumulative hazard function (f(y) = -log(y)) 
and "cloglog" creates a complimentary log-log survival plot  
(f(y) = log(-log(y) along with log scale for the x-axis). 
}
\item{conf.int}{
if \code{TRUE}, confidence bands for the curves are also plotted. 
If set to "only", then only the CI bands are plotted, and the curve 
itself is left off.   
This can be useful for fine control over the colors or line types of a plot. 
}
}
\value{
a list with components x and y, containing the coordinates of the 
last point on each of the curves (but not of the confidence limits). 
This may be useful for labeling. 
}
\section{Side Effects}{
one or more curves are added to the current plot. 
}
\seealso{
\code{\link{lines}},  \code{\link{par}},  \code{\link{plot.survfit}},  \code{\link{survfit}},  \code{\link{survexp}}.   
}
\examples{
# First select by "id" those elements with the greater value of "stop" 
hearta <- by(heart, heart$id, function (x)x[x$stop==max(x$stop),]) 
# Make the resulting list into a data frame by rebinding its 
components 
hearta <- do.call("rbind", hearta) 
# plot the fitted survival curves 
plot(survfit(Surv(stop, event) ~ surgery, data = hearta), lty = 1:2) 
# Add lines denoting the expected survival from National Rate Table 
lines(survexp (stop ~ surgery + ratetable(age = (age + 48) * 365.25, 
      sex = "male", year = year * 365.25), data = hearta), mark = "E",  
      mark.time = 1000, lty = 1:2, cex = 1.5) 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
