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
.AG time2
For interval data only, the ending time of the interval.  Intervals are
assummed to be open on the left and closed on the right, (start, end], and
`event' marks whether an event occured at the end of the interval.
.RT
An object of class 'Surv'.  There are methods for `print', `is.na', and
subscripting survival objects.  To include a survival object inside a
data frame, use the `I()' function.  Surv objects are implimented as
a matrix of 2 or 3 columns.  In time, they will be extended to also include
left censored and interval censored data.
.EX
Surv(aml$time, aml$status)
 [1] 9    13   13+  18   23   28+  31   34   45+  48   161+ 5    5    8    8
 [16] 12   16+  23   27   30   33   43   45
.KW survival
.WR
