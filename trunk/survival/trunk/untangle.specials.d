.BG
.FN untangle.specials
.TL
Help process the 'specials' argument of the terms function.
.DN
Given a terms structure and a desired special name, this return an index
appropriate for subscripting the terms structure and another appropriate for
the data frame.
.CS
untangle.specials(tt, special, order=1)
.RA
.AG tt
a `terms' object.
.AG special
the name of a special function, presumably used in the terms object.
.OA
.AG order
the order of the desired terms.  If set to 2, interactions with the special
function will be included.
.RT
a list with two components:
.AG vars
a vector of variable names, as would be found in the data frame, of the
specials.
.AG terms
a numeric vector, suitable for subscripting the terms structure, that indexes
the terms in the expanded model formula which involve the special.
.EX
# This is code from within the coxph program, m is the data frame and TT
#  the terms structure for the formula.  We wish to exclude first order
#  strata() terms from the final X matrix, and also to extract these terms
#  to a second data frame
TT <- terms(formula, specials='strata')
temp <- untangle.specials(TT, 'strata')
if (length(temp$vars)) {
    X <- model.matrix( TT[-temp$terms], m)
    m.strat <- m[temp$vars]
    }
else X <- model.matrix(TT, m)
.WR
