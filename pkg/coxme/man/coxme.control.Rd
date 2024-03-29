\name{coxme.control}
\Rdversion{1.1}
\alias{coxme.control}
\title{
Auxillary parameters for controlling coxme fits.
}
\description{
Auxillary function which packages the optional parameters of a
\code{coxme} fit as a single list.
}
\usage{
coxme.control(eps = 1e-08, toler.chol = .Machine$double.eps^0.75,
iter.max = 20, inner.iter = Quote(max(4, fit0$iter+1)),
sparse.calc = NULL,
optpar = list(method = "BFGS", control=list(reltol = 1e-5)),
refine.df=4, refine.detail=FALSE, refine.method="control")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eps}{
    convergence criteria for the partial likelihood
}
  \item{toler.chol}{
    tolerance for the underlying Cholesky decomposition.  This is used
    to detect singularity (redundant variables).
}
  \item{iter.max}{
    maximum number of iterations for the final fit
}
  \item{inner.iter}{
    number of iterations for the `inner loop' fits, i.e. when the
    partial likelihood is the objective function of \code{optim}.
    The default is to use one more iteration than the baseline coxph
    model \code{fit0}.  The baseline model contains only the fixed
    effects, and is as part of the setup by the main program.
    The minimum value of 4 applies most often to the case where there
    are no fixed effects.
}
  \item{sparse.calc}{
    choice of method 1 or 2 for a particular portion of the
    calculation.  This can have an effect on run time for problems with
    thousands of random effects.
}
  \item{optpar}{
    parameters passed forward to the \code{optim} routine.
  }
  \item{refine.df}{the degrees of freedom for the t-distribution
    used to draw random samples
    for the refine.n option}
  \item{refine.detail}{this option is mostly for debugging.  If TRUE
    then an extra component \code{refine.detail} will be present in
    the output which contains intermediate variables from the
    iterative refinement calculation.}
  \item{refine.method}{method by which the control calculations are
    done. This is a current research/development question, the option
    will likely disappear at some future date, and users should ignore it.}
}
\details{
 The main flow of \code{coxme} is to use the \code{optim} routine to
 find the best values for the variance parameters.  For any given trial
 value of the variance parameters, an inner loop maximizes the partial
 likelihood to select the regression coefficients beta (fixed) and b
 (random). Within this loop cholesky decomposition is used.  It is
 critical that the convergence criteria of inner loops be less than
 outer ones, thus toler.chol < eps < reltol.
}
\value{
a list of control parameters}
\author{
Terry Therneau
}
\seealso{
\code{\link{coxme}}
}
\keyword{survival}
