.BG
.FN Surv
.TL
Package a survival variable
.CS
Surv(time, event)  or Surv(time, time2, event)
.AG time
for right censored data, this is the follow up time.  For interval data, the
first argument is the starting time for the interval.
.AG event
The status indicator, normally 0=alive, 1=dead.  Other choices are T/F
(TRUE = death) or 1/2 (2=death).
For interval censored data, the status indicator is 1=event at `time',
2= right censored, 3=left censored, 4=interval censored.
.AG time2
For interval censored  or counting process data only, the ending time of the interval.
Intervals are
assummed to be open on the left and closed on the right, (start, end].
For counting process data,
`event' marks whether an event occured at the end of the interval.
.OA
.AG type
one of left, right, counting, or interval.  If this is not specified, the
default is either right or counting, depending on whether the `time2'
argument is absent or present, respectively.
.OA origin
for counting process data, the hazard function origin.  This is most often
used in conjunction with a model containing time dependent strata in order
to align the subjects properly when they cross over from one strata to
another.
.RT
An object of class 'Surv'.  There are methods for `print', `is.na', and
subscripting survival objects.  To include a survival object inside a
data frame, use the `I()' function.  Surv objects are implimented as
a matrix of 2 or 3 columns.
.SH METHOD
In theory it is possible to represent interval censored data without a
third column containing the explicit status.  Exact, right censored,
left censored and interval censored observation would be represented as
intervals of (a,a), (a, infinity), (-infinity,b), and (a,b) respectively;
each specifing the interval within which the event is known to have occured.
Infinity is, of course, impractical in a computer routine.
If the status code is 1, 2 or 3, then the relevant information is assumed to
be contained in `time',  the value in `time2' is ignored, and the second
column of the result will contain a 0.
.pp
At present, all of the methods that handle interval censored data are
parametric models, so the distinction between open and closed intervals
is unimportant.  The distinction is important for counting process data and
the Cox model.
.EX
Surv(aml$time, aml$status)
 [1] 9    13   13+  18   23   28+  31   34   45+  48   161+ 5    5    8    8
 [16] 12   16+  23   27   30   33   43   45
.KW survival
.WR
