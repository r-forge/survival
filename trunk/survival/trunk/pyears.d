.BG
.FN pyears
.TL
Person Years
.DN
This function computes the person-years of follow-up time contributed by a
cohort of subjects, stratifed into subgroups.
It also computes the number of subjects who contribute to each cell of the
output table, and optionally the number of events and/or expected number of
events in each cell.
.CS
pyears(formula, data, weights, subset, na.action, ratetable=survexp.us,
scale=365.25, model=F, x=F, y=F)
.RA
.AG formula
a formula object.  The response variable will be a vector of follow-up times
for each subject, or a Surv object containing the follow-up time and an
event indicator.
The predictors consist of optional grouping variables
separated by + operators (exactly as in `survfit'), time-dependent grouping
variables such as age (specified with `tcut'), and optionally a
`ratetable()' term.  This latter matches each subject to his/her expected
cohort.
.OA
.AG data
a data frame in which to interpret the variables named in the formula.
.AG weights
case weights.
.AG subset
expression stating that only a subset of the rows should be used.
.AG na.action
missing data filter function.
.AG ratetable
a table of event rates, such as survexp.uswhite.
.AG scale
a scaling for the results.  As most rate tables are in units/day, the
default value of 365.25 causes the output to be reported in years.
.AG "model, x, y"
If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result.
.RT
a list with components
.RC pyears
an array containing the person-years of exposure. (Or other units, depending
on the rate table and the scale).
.RC n
an array containing the number of subjects who contribute time to each cell
of the pyears array.
.RC event
an array containing the observed number of events.  This will be present only
if the resonse variable is a Surv object.
.RC expeced
an array containing the expected number of events.  This will be present only
if there was a ratetable term.
.RC offtable
the number of person-years of exposure in the cohort that was not part of
any cell in the pyears array.  This is often useful as an error check; if
there is a mismatch of units between two variables, nearly all the person
years may be off table.
.RC call
an image of the call to the function.
.RC na.action
the na.action attribute contributed by an na.action routine, if any.
.DT
Because pyears may have several time variables, it is necessary that all
of them be in the same units.  For instance in the call
.Cs
 py <- pyears(futime ~ rx + ratetable(age=age, sex=sex, year=entry.dt))
.Ce
with a ratetable whose natural unit is days, it is important that futime,
age and entry.dt all be in days.  Given the wide range of possible inputs,
it is difficult for the routine to do sanity checks of this aspect.
.pp
A special function `tcut' is needed to specify time-dependent cutpoints.
For instance, assume that age is in years, and that the desired final
arrays have as one of their margins the age groups 0-2, 2-10, 10-25, and 25+.
A subject who enters the study at age 4 and remains under observation for
10 years will contribute follow-up time to both the 2-10 and 10-25
subsets.  If `cut(age, c(0,2,10,25,100))' were used in the formula, the
subject would be classifed according to his starting age only.  The tcut
function has the same arguments as cut, but produces a different output
object which allows the pyears function to correctly track the subject.
.pp
The results of pyears() are normally used as input to further calculations.
The print routine, therefore, is designed to give only a summary of the
table.
The example below is from a study of hip fracture rates from 1930 - 1990
in Olmstead County, Minnesota.  Survival post hip fracture has increased over
that time, but so has the survival of elderly subjects in the population at
large.  A model of relative survival helps to clarify what has happened:
Poisson regression is used, but replacing exposure time with expected
exposure (for an age and sex matched control).
Death rates change with age, of course, so the result is carved into
1 year increments of time.  Males and females were done separately.
.EX
attach(hips)
temp1 <- tcut(dt.fracture, seq(from=mdy.date(1,1,30), by=365.25, length=61)
temp2 <- tcut(age*365.5,   365.25*(0:105))   #max age was > 100!
pfit.m<- pyears(Surv(futime, status) ~ temp1 + temp2 +
			 ratetable(age=age*365.25, year=dt.fracture, sex=1),
	       subset=(sex==1),
	       ratetable=survexp.minnwhite)
# now, convert the arrays into a data frame
tdata <- data.frame( age  = (0:105)[col(pfit.m$pyears)],
		     yr   = (1930:1990)[row(pfit.m$pyears)],
		       y  = c(pfit.m$event),
		     time = c(pfit.m$expect))
# fit the gam model
gfit.m <- gam(y ~ s(age) + s(yr) + offset(log(time)), family=poisson,
			data= tdata)
plot(gfit.m, se=T)
.SA
ratetable, survexp, Surv
.KW survival
.WR
