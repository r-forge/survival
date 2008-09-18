\name{survfit}
\alias{survfit}
\title{
  Compute a Survival Curve for Censored Data 
}
\description{
Computes an estimate of a survival curve for censored data 
using either the Kaplan-Meier or the Fleming-Harrington method 
or computes the predicted survivor function for a Cox proportional 
hazards model. 
}
\usage{
survfit(formula, data=sys.parent(), weights, subset, na.action,  
        newdata, individual=F, conf.int=.95, se.fit=T,  
        type=&lt;&lt;see below&gt;&gt;, error=&lt;&lt;see below&gt;&gt;, 
        conf.type="log", conf.lower="usual", start.time)
}
\arguments{
\item{formula}{
a formula object or a \code{coxph} object. 
If a formula object is supplied it must have a 
\code{Surv} object as the  
response on the left of the ~ operator and, if desired, terms  
separated by + operators on the right. 
One of the terms may be a \code{strata} object.
}
\item{data}{
a data frame in which to interpret the variables named in the formula, 
or in the \code{subset} and the \code{weights} argument. 
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
\item{na.action}{
a missing-data filter function, applied to the model frame, after any 
\code{subset} argument has been used. 
Default is \code{options()\$na.action}. 
}
\item{newdata}{
a data frame with the same variable names as those that appear 
in the \code{coxph} formula.  
Only applicable when \code{formula} is a 
\code{coxph} object. 
The curve(s) produced will be representative of a cohort who's 
covariates correspond to the values in \code{newdata}. 
Default is the mean of the covariates used in the 
\code{coxph} fit. 
}
\item{individual}{
a logical value, if \code{TRUE}, 
the data frame represents different 
time epochs for only one individual.  
If \code{FALSE}, multiple rows indicate 
multiple individuals.  
If \code{TRUE}, only one curve 
will be produced, 
if \code{FALSE}, there will be one curve per row in 
\code{newdata}.  
Only applicable when \code{formula} is 
a \code{coxph} object. 
}
\item{conf.int}{
the level for a two-sided confidence interval on the survival curve(s). 
Default is 0.95. 
}
\item{se.fit}{
a logical value indicating whether standard errors should be 
computed.  Default is \code{TRUE}. 
}
\item{type}{
a character string specifying the type of survival curve. 
Possible values are \code{"kaplan-meier"}, 
\code{"fleming-harrington"} or 
\code{"fh2"} 
if a formula is given 
and \code{"aalen"} or 
\code{"kaplan-meier"} 
if the first argument is a \code{coxph} object, 
(only the first two characters are necessary). 
The default is \code{"aalen"} when 
a \code{coxph} object is given, 
and it is \code{"kaplan-meier"} otherwise.  
Earlier versions of \code{survfit} 
used \code{type="tsiatis"} to 
get the \code{"aalen"} estimator.  For backward 
compatibility, this is still allowed. 
}
\item{error}{
a character string specifying the error.  Possible values are 
\code{"greenwood"} for the Greenwood formula or 
\code{"tsiatis"} for the Tsiatis formula, 
(only the first character is 
necessary).  
The default is \code{"tsiatis"} when 
a \code{coxph} object is 
given, and it is \code{"greenwood"} otherwise. 
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
\item{start.time}{
numeric value specifying a time to start calculating survival information. 
}
}
\value{
an object of class \code{"survfit"}.  
See \code{survfit.object} for 
details. Methods defined for survfit objects are  
\code{print}, \code{plot}, 
\code{lines}, and \code{points}. 
}
\details{
The estimates used are the Kalbfleisch-Prentice 
(Kalbfleisch and Prentice, 1980, p.86) and the Tsiatis/Link/Breslow, 
which reduce to the Kaplan-Meier and Fleming-Harrington estimates, 
respectively, when the weights are unity.  
When the parent data is a Cox model, there is an extra
term in the variance of the curve, due to the variance of the coefficients and 
hence variance in the computed weights.  Details of the Aalen 
estimator and its variance are found in Tsiatis (1981). 

The Greenwood formula for the variance is a sum of terms 
d/(n*(n-m)), where d is the number of deaths at a given time point, n 
is the sum of weights for all individuals still at risk at that time, and 
m is the sum of weights for the deaths at that time.  The 
justification is based on a binomial argument when weights are all 
equal to one; extension to the weighted case is ad hoc.  Tsiatis 
(1981) proposes a sum of terms d/(n*n), based on a counting process 
argument which includes the weighted case. 

For the Fleming-Harrington estimate, two different methods for 
handling ties have been implemented. 
If there were 3 deaths out of 10 at risk, then the original estimate
of Nelson increments the hazard by 3/10, while the modification of
Fleming and Harrington increments 
it by 1/10 + 1/9 + 1/8.  For curves created 
after a Cox model these correspond to the Breslow and Efron estimates, 
respectively, and the proper choice is made automatically. 
The \code{fh2} method will 
give results closer to the Kaplan-Meier. 

When the data set includes left censored or interval censored data (or both),
then the EM approach of Turnbull is used to compute the overall curve.
When the baseline method is the Kaplan-Meier, this is known to converge to
the maximum likelihood estimate.

Based on the work of Link (1984), the log transform is expected to produce 
the most accurate confidence intervals.  If there is heavy censoring, then 
based on the work of Dorey and Korn (1987) the modified estimate will give 
a more reliable confidence band for the tails of the curve. 
}
\section{References}{
Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for confidence 
intervals for survival probabilities.  \emph{Statistics in Medicine} 
\bold{6}, 679-87. 
Fleming, T. H. and Harrington, D. P. (1984).  Nonparametric estimation of the 
survival distribution in censored data.  \emph{Comm. in Statistics}  
\bold{13}, 2469-86. 
Kablfleisch, J. D. and Prentice, R. L. (1980).  
\emph{The Statistical Analysis of Failure Time Data.}
New York:Wiley. 
Link, C. L. (1984). Confidence intervals for the survival 
function using Cox's proportional hazards model with  
covariates.  \emph{Biometrics}  
\bold{40}, 601-610.
Turnbull, B. W. (1974).  Nonparametric estimation of a survivorship
function with doubly censored data. \emph{J Am Stat Assoc},
\bold{69}, 169-173. 
Tsiatis, A. (1981). A large sample study of the estimate 
for the integrated hazard function in Cox's regression 
model for survival data. \emph{Annals of Statistics}  
\bold{9}, 93-108. 
}
\seealso{
\code{\link{print}},  
\code{\link{plot}},  
\code{\link{lines}},   
\code{\link{coxph}},  
\code{\link{Surv}},  
\code{\link{strata}}.   
}
\examples{
#fit a Kaplan-Meier and plot it 
fit <- survfit(Surv(time, status) ~ group, data = leukemia) 
plot(fit, lty = 2:3) 
legend(100, .8, c("Maintained", "Nonmaintained"), lty = 2:3) 
#fit a Cox proportional hazards model and plot the  
#predicted survival curve at the average predictor 
fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian) 
plot(survfit(fit), xlab = "Survival in Days") 

# Here is the data set from Turnbull
#  There are no interval censored subjects, only left-censored (status=3),
#  right-censored (status 0) and observed events (status 1)
tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
                    status=rep(c(1,0,2),4),
                    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))
fit  <- survfit(Surv(time, time, status, type='interval') ~1, 
              data=tdata, weight=n)
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
