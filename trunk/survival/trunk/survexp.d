.BG
.FN survexp
.TL
Compute expected survival
.DN
Computes the expected survival post study entry for a cohort of subjects,
based on their age and sex composition.
.CS
survexp(entry, birth, sex,
data, subset, na.action,
times=round(182.6 * 0:8),
type=c("mean", "individual", "matrix"), expected=survexp.uswhite, interp=F)
.RA
.AG entry
date of entry to the study.  Must be a date value as returned by mdy.date.
.AG birth
date of birth, a date value as returned by mdy.date.  Optionally, this may
contain the patient's age, in which case the date of birth is inferred by
subtracting round(365.25*age) from the entry date.
.AG sex
this indexes the second dimension of the expected array. For the example
provided-- US white-- the values are 1=male, 2=female.
.OA
.AG data
frame in which to find the variables.
.AG subset, na.action
the subset of observations to be used in the calculation, and the function
to be used to handle any NAs in the data.
.AG times
vector of time points at which an expected survival is desired.  If survexp
is used for plotting purposes, then this can be a fairly coarse grid, since
the expected survival function is very smooth.
.AG type
if 'mean', then the average expected survival for the cohort is returned.
The result will often be added to an existing Kaplan-Meier plot.
If type=='individual', then the times vector must be of the same length as
entry and birth date, and the return value contains the expected survival
of each individual at his/her respective time.
This form is most often used with the one sample log-rank test, see survdiff.
If type=='matrix', then one column is returned per subject.
.AG expected
a 2 or 3 dimensional array, age by sex by calendar year, containing the
expected survival data.
The second dimension can be of any size >1.  For instance, a table containing
both sex and 4 catorizations of race would have second dimension 8.
The `sex' argument is understood to be an index into this dimension of the
table.
.AG interp
if T, then the `expected' matrix is filled out to contain one entry per year,
by linear interpolation.
.RT
a structure with components
.AG  time
vector of time values
.AG surv
if mean=T, a vector of survival probablities for the cohort as a whole.  If
mean=F and length(times) = number of subjects, then y is a vector of
survival probabilities for each subject at that subject's time value.
Otherwise y will be a matrix containing each subject's expected survival at
each time in the times vector.
.AG n
number of subjects, after removal of missing values.  An attribute 'omit'
contains the indices of the removed observations, if any.
.SH METHOD
The expected matrix is constructed from national life tables; a
sample of references is below.
We make use of the life table `q' for each year of age, which is the
probability that a subject who reaches that age will die before reaching
age+1.  For example, in the 1980 table for white males, the row for age
62 is the probability that someone who becomes 62 in 1980 will die before his
63rd birthday.
.pp
To make interpolation easier, the expected table used by the program contains
the per day hazard, i.e., -log(1-q) /365.25.  For a 62 year old male in 1980
q=.0212 and h= .0004868.  The cumulative hazard over an
interval is then just #days * hazard.
For each subject, a cumulative hazard
function is constructed which is piecewise linear in intervals starting at
their birth date.  The values at the desired times from entry to the study
are then plucked off of this curve, and the survival at those times
is computed as exp(-cum hazard).  The overall survival is the mean of the
individual patient survivals.
.pp
The dimnames of the hazard array `expected' control which entries of the
array are used.  The first name vector contains the starting age for each
age interval.  The first interval must start at 0, but intervals need not be
sequential integers as used in the default survexp.uswhite data set, e.g.,
a study of infant mortality would likely wish to divide the first year more
finely, or a non-US population might have data only on 5 year age groups.
(The internal calculations are done in days, so there is a lower limit to the
granularity.)  The third subscript contains the starting year.  For instance,
if 60-65 were an age interval, then the calendar year in which the subject
turns 60 will be used in choosing the appropriate hazard rate from the array.
This hazard rate is applied for the entire age interval.
.pp
The last age and the last calendar year of the array are both considered to
be open ended, i.e., they apply to all later ages and calendar years,
respectively.  The first calendar year is applied for all earlier years
(in defiance of it's labeled starting date).
.SH REFERENCES
National Office of Vital Statistics: Life Tables for the Geographic Divisions
of the United States: 1959-61.  Vol. 1, No. 3.  Public Health Service,
Washington.  U.S. Government Printing Office, May 1965.
.SA
survfit, survdiff
.EX
xx <-  survfit(last.dt - entry.dt, last.status)
plot(xx)
maxt <- max(xx$time)
lines(survexp(entry.dt, birth.dt, sex, times=seq(0,maxt,length=10), lty=2)
.KW survival
.WR
