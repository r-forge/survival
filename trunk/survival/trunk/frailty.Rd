\name{frailty}
\alias{frailty}
\title{
  Random effects terms
}
\description{
The frailty function allows one to add a simple random effects term to
a Cox or survreg model.
}
\usage{
frailty(x, distribution="gamma", sparse, theta, df, method, eps)
}
\arguments{
\item{x}{
the variable to be entered as a random effect.  
It is always treated as a factor.
}
\item{distribution}{
either the \code{gamma}, 
\code{gaussian} or \code{t}
distribution may be specified.
The routines \code{frailty.gamma},
\code{frailty.gaussian} and 
\code{frailty.t} do the actual work.
}
\item{sparse}{
cutoff for using a sparse coding of the data matrix.  
If the total number of levels of \code{x} is larger
than this value, then sparse coding is used.
The correct cutoff is still a matter of exploration: if the number of
levels is very large (thousands) then the non-sparse calculation may not be
feasable in terms of both memory and compute time.  
Likewise, the accuracy of the sparse approximation appears to be related to
the maximum proportion of subjects in any one class, being best when no
one class has a large membership.
}
\item{theta}{
if specified, this fixes the variance of the random effect.
If not, the variance is a parameter, and a best solution is sought.
Specifying this implies \code{method='fixed'}.
}
\item{df}{
if specified, this fixes the degrees of freedom for the random effect.
Specifying this implies \code{method='df'}.
Only one of \code{theta} or 
\code{df} should be specified.
}
\item{method}{
the method used to select a solution for theta, the variance of the
random effect.  
The \code{fixed} corresponds to a user-specified
value, and no iteration is done.
The \code{df} selects the variance such that the
degrees of freedom for the random effect matches a user specified value.
The \code{aic} method seeks to 
maximize Akiake's information criteria 
2*(partial likelihood - df).
The \code{ml} and \code{reml}
methods are specific to Cox models with gamma and gaussian random effects,
respectively.
Please see further discussion below.
}
\item{tdf}{
the degrees of freedom for the t-distribution.
}
\item{eps}{
convergence critera for the iteration on theta.
}
}
\value{
this function is used in the model statment of either
\code{coxph} or \code{survreg}.
It's results are used internally.
}
\details{
The \code{frailty} plugs into the general penalized
modeling framework provided by the \code{coxph} 
and \code{survreg} routines.  
This framework deals with likelihood, penalties, and degrees of freedom;
these aspects work well with either parent routine.

Therneau, Grambsch, and Pankratz show how maximum likelihood estimation for
the Cox model with a gamma frailty can be accomplished using a general
penalized routine, and Ripatti and Palmgren work through a similar argument
for the Cox model with a gaussian frailty.  Both of these are specific to
the Cox model.  
The first in particular has a flavor of `fooling' the penalized code into
doing the needed calculation.
Use of gamma/ml or gaussian/reml with 
\code{survreg} does not lead to valid results. 

The extensible structure of the penalized methods is such that the penalty
function, such as \code{frailty} or
\code{pspine}, is completely separate from the modeling
routine.  The strength of this is that a user can plug in any penalization
routine they choose.  A weakness is that it is very difficult for the
modeling routine to know whether a sensible penalty routine has been
supplied.  
}
\section{References}{
S Ripatti and J Palmgren, Estimation of multivariate frailty 
models using penalized partial likelihood,
Biometrics, 56:1016-1022, 2000.

T Therneau, P Grambsch and VS Pankratz,
Penalized survival models and frailty,
J Computational and Graphical Statistics, 12:156-175, 2003.
}
\seealso{
coxph, survreg
}
\examples{
# Random institutional effect
> coxph(Surv(time, status) ~ age + frailty(inst, df=4), lung)

                        coef se(coef)     se2 Chisq   DF     p 
                  age 0.0194 0.00933  0.00925 4.31  1.00 0.038
frailty(inst, df = 4)                         3.34  3.99 0.500

Iterations: 3 outer, 8 Newton-Raphson
     Variance of random effect= 0.0379   I-likelihood = -743.6 
Degrees of freedom for terms= 1 4 
Likelihood ratio test=9.96  on 4.97 df, p=0.075
  n=227 (1 observations deleted due to missing values)
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
