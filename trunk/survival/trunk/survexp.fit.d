.BG
.FN survexp.fit
.TL
Compute expected survival
.CS
survexp.fit(x, y, times, ratetable, rfac)
.RA
.AG x
The first column contains the group, an integer value that divides the
subjects into subsets.  Remaining columns must match the dimensions of
the ratetable, in the correct order.
.AG y
the follow up time for each subject.
.AG times
the vector of times at which a result will be computed.
.AG death
indictes whether or not `y' includes death times.
.AG ratetable
a rate table, such as survexp.uswhite.
.RT
A list containing the number of subjects, a sum of weights,
and the weighted sum of expected survivals for each time interval.
If there are multiple groups, these will be
matrices with one column per group.
.DT
For any time interval (times[i], times[i+1]), n will be the number of subjects
who are at risk for some or all of the interval, i.e., all those for which
y > times[i].  Let h[j], j= 1 to n, be the hazard experienced
during the interval by each of these subjects,
and w[j] be the weight, where w[j] = Pr(survival to the beginning of the
inteval).
The survival for the interval is defined as the weighted mean
of exp(-h).  If n>1, this average makes sense only if all
subjects are at risk for the entire interval, which means that `times'
should be a superset of `y'.
(If not, then an estimate that accounted for drop-outs within the interval
should be used; the routine is not yet that smart.)
.SH WARNING
Most users will call the higher level routine `survexp'.
Consequently, this function has very few error checks on its input arguments.
.SA
survexp, survexp.uswhite
.KW survival
.WR
