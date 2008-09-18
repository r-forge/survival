\name{survexp}
\alias{survexp}
\title{
  Compute Expected Survival 
}
\description{
Returns either the expected survival of a cohort of subjects, or the 
individual expected survival for each subject. 
}
\usage{
survexp(formula, data, weights, subset, na.action, times, cohort=T, 
        conditional=F, ratetable=survexp.us, scale=1, npoints, 
        se.fit=&lt;&lt;see below&gt;&gt;, model=F, x=F, y=F)
}
\arguments{
\item{formula}{
formula object.  The response variable is a vector of follow-up times 
and is optional.  The predictors consist of optional grouping variables 
separated by the \code{+} operator (as in \code{survfit}), along with a \code{ratetable}  
term.  The \code{ratetable} term matches each subject to his/her expected cohort. 
}
\item{data}{
data frame in which to interpret the variables named in 
the \code{formula}, \code{subset} and \code{weights} arguments. 
}
\item{weights}{
case weights. 
}
\item{subset}{
expression indicating a subset of the rows of \code{data} to be used in the fit. 
}
\item{na.action}{
function to filter missing data. This is applied to the model frame after  
\code{subset} has been applied.  Default is \code{options()\$na.action}. A possible 
value for \code{na.action} is \code{na.omit}, which deletes observations that contain 
one or more missing values. 
}
\item{times}{
vector of follow-up times at which the resulting survival curve is  
evaluated.  If absent, the result will be reported for each unique  
value of the vector of follow-up times supplied in \code{formula}. 
}
\item{cohort}{
logical value: if \code{FALSE}, each subject is treated as a subgroup of size 1. 
The default is \code{TRUE}. 
}
\item{conditional}{
logical value: if \code{TRUE}, the follow-up times supplied in \code{formula} 
are death times and conditional expected survival is computed. 
If \code{FALSE}, the follow-up times are potential censoring times.  
If follow-up times are missing in \code{formula}, this argument is ignored.   
}
\item{ratetable}{
table of event rates, such as \code{survexp.uswhite}, or a fitted Cox model. 
}
\item{scale}{
numeric value to scale the results.  If \code{ratetable} is in units/day, 
\code{scale = 365.25} causes the output to be reported in years. 
}
\item{npoints}{
number of points at which to calculate intermediate results, evenly spaced  
over the range of the follow-up times.  The usual (exact) calculation is done  
at each unique follow-up time. For very large data sets specifying \code{npoints}  
can reduce the amount of memory and computation required. 
For a prediction from a Cox model \code{npoints} is ignored. 
}
\item{se.fit}{
logical value indicating whether to compute the standard error of the  
predicted survival.  The default is to compute standard errors whenever  
possible, which, at this time, is only for the Ederer method and a Cox  
model as the rate table. 
}
\item{model}{
logical value: if \code{TRUE}, the model frame is included as component 
\code{model} in the return object.  
}
\item{x}{
logical value: if \code{TRUE}, the model matrix is returned as component 
\code{x} in the return object. 
}
\item{y}{
logical value: if \code{TRUE}, the response in \code{formula} is returned as component  
\code{y} in the return object.  
}
}
\value{
if \code{cohort=TRUE} an object of class \code{survexp}, 
otherwise a vector of per-subject expected survival values. 
The former contains the number of subjects at risk 
and the expected survival for the cohort at each requested time. 
}
\details{
Individual expected survival is usually used in models or testing, to 
"correct" for the age and sex composition of a group of subjects. 
For instance, assume that birth date, entry date into the study, 
sex and actual survival time are all known for a group of subjects. 
The \code{survexp.uswhite} population tables contain expected death rates 
based on calendar year, sex and age. 
Then 
<s-expression>
haz <- -log(survexp(death.time ~ ratetable(sex=sex, year=entry.dt, 
            age=(birth.dt-entry.dt)), cohort=F)) 
</s-expression>
gives for each subject the total hazard experienced up to their observed 
death time or censoring time. 
This probability can be used as a rescaled time value in models: 
<s-expression>
glm(status ~ 1 + offset(log(haz)), family=poisson) 
glm(status ~ x + offset(log(haz)), family=poisson) 
</s-expression>
In the first model, a test for intercept=0 is the one sample log-rank 
test of whether the observed group of subjects has equivalent survival to 
the baseline population.  The second model tests for an effect of variable 
\code{x} after adjustment for age and sex. 

Cohort survival is used to produce an overall survival curve.  This is then 
added to the Kaplan-Meier plot of the study group for visual comparison 
between these subjects and the population at large.  There are three common 
methods of computing cohort survival. 
In the "exact method" of Ederer the cohort is not censored; this corresponds 
to having no response variable in the formula.  Hakulinen recommends censoring 
the cohort at the anticipated censoring time of each patient, and Verheul 
recommends censoring the cohort at the actual observation time of each 
patient. 
The last of these is the conditional method. 
These are obtained by using the respective time values as the 
follow-up time or response in the formula. 
}
\section{References}{
Berry, G. (1983). The analysis of mortality by the subject-years method. 
\emph{Biometrics}
\bold{39},173-84. 

Ederer, F., Axtell, L. and Cutler, S. (1961). 
The relative survival rate: a statistical methodology. 
\emph{Natl Cancer Inst Monogr}
\bold{6},101-21. 

Hakulinen, T. (1982). 
Cancer survival corrected for heterogeneity in patient withdrawal. 
\emph{Biometrics}
\bold{38}, 933. 

Verheul, H., Dekker, E., Bossuyt, P., Moulijn, A. and Dunning, A. (1993). 
Background mortality in clinical survival studies. 
\emph{Lancet}
\bold{341}, 872-875. 
}
\seealso{
\code{\link{by}},  \code{\link{do.call}},  \code{\link{pyears}},  \code{\link{survexp.us}},  \code{\link{survexp.fit}},  \code{\link{survfit}}.   
}
\examples{
# Create new data frame with the largest stop value for each patient 
hearta <- by(heart, IND=heart$id, FUN=function(x)x[x$stop==max(x$stop),]) 
hearta <- do.call("rbind", hearta) 
# Estimate of conditional survival  
survexp(stop ~ ratetable(sex="male", year=year*365.25,  
          age=(age+48)*365.25), conditional=T, data=hearta) 
# Estimate of conditional survival stratified by prior surgery 
survexp(stop ~ surgery + ratetable(sex="male", year=year*365.25,  
	age=(age+48)*365.25), conditional=T, data=hearta) 
rm(hearta) 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
