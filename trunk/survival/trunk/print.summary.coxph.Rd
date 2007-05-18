\name{print.summary.coxph}
\alias{print.summary.coxph}
\title{
  Print method for summary.coxph objects
}
\description{
Produces a printed summary of a fitted coxph model 
}
\usage{
summary.coxph(x, coef=T, conf.int=T, digits=max(options()$digits - 4, 3))
}
\arguments{
\item{x}{
the result of a call to \code{summary.coxph} 
}
\item{coef}{
logical flag: print the table of coefficients? 
}
\item{conf.int}{
logical flag: print the table of confidence intervals?
}
