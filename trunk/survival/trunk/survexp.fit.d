.BG
.FN survexp.fit
.TL
Compute expected survival
.CS
survexp.fit(x, y, times, conditional, ratetable)
.RA
.AG x
The first column contains the group, an integer value that divides the
subjects into subsets.  Remaining columns must match the dimensions of
the ratetable, in the correct order.
.AG y
the follow up time for each subject.
.AG times
the vector of times at which a result will be computed.
.AG conditional
if T compute the conditional survival, if false compute cohort survival.
.AG ratetable
a rate table, such as survexp.uswhite.
.RT
A list containing the number of subjects and the expected survival(s)
at each time point.
If there are multiple groups, these will be
matrices with one column per group.
.DT
For conditional survival y must be the time of last follow-up or death for
each subject.  For cohort survival it must be the potential censoring time for
each subject, ignoring death.
.pp
For an exact estimate `times' should be a superset of `y', so that each
subject at risk is at risk for the entire sub-interval of time.
For a large data set, however, this can use an inordinate amount of
storage and/or compute time.  If the `times' spacing is more coarse than
this, an actuarial approximation is used which should, however, be extremely
accurate as long as all of the returned values are > .99.
.SH WARNING
Most users will call the higher level routine `survexp'.
Consequently, this function has very few error checks on its input arguments.
.SA
survexp, survexp.us
.KW survival
.WR
