.BG
.FN is.ratetable
.TL
Verify that an object is of class `ratetable'.
.DN
The function verifies not only the `class' attribute, but the
structure of the object.
.CS
is.ratetable(x, verbose=F)
.RA
.AG x
the object to be verified.
.OA
.AG verbose
if TRUE and the object is not a ratetable, then return a character string
describing the way(s) in which `x' fails to be a proper ratetable object.
.RT
returns TRUE if x is a ratetable, and FALSE or a description if it is not.
.DT
Rate tables are used by the `pyears' and `survexp' functions, and normally
contain death rates for some population, categorized by age, sex, or other
variables.  They have a fairly rigid structure, and the `verbose' option
can help in creating a new rate table.
.SA
pyears, survexp
.KW survival
.WR
