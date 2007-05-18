\name{untangle.specials}
\alias{untangle.specials}
\title{
  Process the <s-expression>specials</s-expression> Argument of the Terms Function 
}
\description{
Given a terms structure and a desired special name, this returns an index 
appropriate for subscripting the terms structure and another appropriate for 
the data frame. 
This function is used by the survival functions. 
}
\usage{
untangle.specials(tt, special, order=1)
}
\arguments{
\item{tt}{
a \code{terms} object. 
}
\item{special}{
the name of a special function, presumably used in the terms object. 
}
\item{order}{
the order of the desired terms. 
If set to 2, interactions with the special function will be included. 
}
}
\value{
a list with two components: 

\item{vars}{
a vector of variable names, as would be found in the data frame, of the 
specials. 
}
\item{terms}{
a numeric vector, suitable for subscripting the terms structure, that indexes 
the terms in the expanded model formula which involve the special. 
}
}
\examples{
# This is code from within the coxph program, m is the data frame and TT 
# the terms structure for the formula.  We wish to exclude first order 
# strata() terms from the final X matrix, and also to extract these terms 
# to a second data frame 
TT <- terms(formula, specials='strata') 
temp <- untangle.specials(TT, 'strata') 
if (length(temp$vars)) { 
    X <- model.matrix( TT[-temp$terms], m) 
    m.strat <- m[temp$vars] 
    } 
else X <- model.matrix(TT, m) 
}
\keyword{survival4}
% docclass is function
% Converted by Sd2Rd version 37351.
