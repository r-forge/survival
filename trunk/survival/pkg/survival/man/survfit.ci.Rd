\name{survfit.ci}
\alias{survfit.ci}
\title{
  Compute the cumulative incidence function
}
\description{
In competing risk models the cumulative incidence function is a
compelling alternative to the usual Kaplan-Meier.
This function also compute the current incidence functions for
multiple-state models.
}
\usage{
survfit.ci(formula=formula(data), data=sys.parent(), weights, 
     subset, na.action, absorb, 
     type=c("kaplan-meier", "fleming-harrington", "fh2"), 
     conf.int=0.95, se.fit=T, conf.type=c("log", "log-log", "plain", "none"), 
     conf.lower=c("usual", "peto", "modified"), model=F, x=F, y=F)
}
\arguments{
\item{formula}{
a formula object or a \code{coxph} object. 
If a formula object is supplied it must have a 
\code{Surv} object as the  
response on the left of the ~ operator and, if desired, terms  
separated by + operators on the right. 
The right hand side must contain a \code{strata}
expression.
}
\item{data}{
a data frame in which to interpret the variables named in the formula, 
or in the \code{subset} and the 
\code{weights} argument. 
}
\item{weights}{
case weights, they must be nonnegative and it is strongly recommended that  
they be strictly positive, since zero weights are ambiguous, compared 
to use of the \code{subset} argument. 
}
\item{subset}{
expression saying that only a subset of the rows of the data 
should be used in the fit. 
}
\item{na.action }{
a missing-data filter function, applied to the model frame, after any 
\code{subset} argument has been used. 
Default is \code{options()\$na.action}. 
}
\item{absorb}{
for current incidence calculations, the list of absorbing states
}
\item{type}{
a character string specifying the type of survival curve. 
Possible values are \code{"kaplan-meier"}, 
\code{"fleming-harrington"} or 
\code{"fh2"}. 
}
\item{conf.type}{
a character string specifying the confidence interval type.  Possible 
values are: \code{"none"} for no confidence intervals, 
\code{"plain"} for the
standard intervals, "curve +- k *se(curve)", where k is determined 
from \code{conf.int},  
\code{"log"} for intervals based on the 
cumulative hazard 
or log(survival) (the default), 
and \code{"log-log"} for intervals based on 
the log hazard or log(-log(survival)).  The last type will never 
extend past 0 or 1.  Only 
enough of the string to uniquely identify it is necessary. 
}
\item{conf.lower}{
a character string to specify modified lower limits to the curve, the 
upper limit remains unchanged.  
Possible values are \code{"usual"} (unmodified), 
\code{"peto"}, 
and \code{"modified"}.  T
he modified lower limit 
is based on an "effective n" argument.  The confidence 
bands will agree with the usual calculation at each death time, but unlike 
the usual bands the confidence interval becomes wider at each censored 
observation.  The extra width is obtained by multiplying the usual 
variance by a factor m/n, where n is the number currently at risk and 
m is the number at risk at the last death time.  (The bands thus agree 
with the un-modified bands at each death time.) 
This is especially useful for survival curves with a long flat tail. 

The Peto lower limit is based on the same "effective n" argument as the 
modified limit, but also replaces the usual Greenwood variance term with 
a simple approximation.  It is known to be conservative. 
}
\item{conf.int }{
the level of the confidence interval
}
}
\value{
a cumulative incidence estimate object.
This structure is identical to that for a 
\code{survfit} object, with the exception that
the \code{surv} component will be a matrix
with multiple columns, one for each competing cause.
}
\details{
The cumulative incidence curve is an alternative to the Kaplan-Meier
for competing risks data.
For instance, in patients with MGUS, conversion to an overt plasma
cell malignancy occurs at about 1\% per year.  
A Kaplan-Meier estimate, treating death due to other causes as censored,
gives a 20 year cumulate rate of 33\% for the 241 early patients of Kyle.  
This estimates the incidence of conversion, if other causes of death were
removed.

  
The CI estimate, on the other hand, estimates the total number of
conversions that will actually occur.  Because the population is older,
this is much smaller than the KM, 22\% at 20 years for Kyle's data.
If there were no censoring, then CI(t) would simply be the total number
of patients with progression by time t, divided by the sample size n.

  
On the left of the equation, the \code{status}
variable indicates subjects who were censored or had an event of any
type, as in the usual survival functons.
On the right hand side, \code{strata()} 
indicates the type of event.  (The value of the strata is ignored for
subjects who are censored, but do not set it to NA, as that would cause
all the censored observations to be deleted!).
Other variables on the right hand side define subgroups, as they
do for \code{survfit}.

For multi-state models, a given subject will have multiple (start, stop]
records covering their time under observation, for each of these
the \code{strata} identifies the state to which
they transfered at the end of the interval.
A \code{cluster(id)} term names an identifier variable
that links the observations for a single individual together.
In this case, \code{survfit.ci} estimates the
total fraction of subjects in each state over time.
The standard error in this case is an approximation whose accuracy has not
yet been widely explored.
}
\section{References}{
Kablfleisch, J. D. and Prentice, R. L. (1980).  
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley.

Kyle, R. A. (1997).
Moncolonal gammopathy of undetermined significance and solitary
      plasmacytoma. Implications for progression to overt multiple myeloma\},
\emph{Hematology/Oncology Clinics N. Amer.}
\bold{11}, 71-87.
}
\seealso{
\code{\link{survfit}}  
}
\examples{
# Separate curves for males and females, with endpoints of death
#  or progression to a plasma cell malignancy
endpoint <- ifelse(is.na(mgus$sdtime), 'death', 'progression')
fit <- survfit.ci(Surv(time, status) ~ sex + strata(endpoint), mgus)
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
