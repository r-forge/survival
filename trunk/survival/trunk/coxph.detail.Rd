\name{coxph.detail}
\alias{coxph.detail}
\title{
  Details of a Cox Model Fit 
}
\description{
Returns the individual contributions to the first and second derivative 
matrix, at each unique event time. 
}
\usage{
coxph.detail(fit)
}
\arguments{
\item{fit}{
a Cox model object, i.e., the result of \code{coxph}. 
}
\item{riskmat}{
should the matrix of "at risk" indictators be returned
}
