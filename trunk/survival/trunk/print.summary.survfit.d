.\" @(#)prt_sm_srvft.d version 3.1 created 4/5/96
.\" @(#)Copyright (c), 1987, 1996 StatSci, Inc.  All rights reserved.
.BG
.FN print.summary.survfit
.TL
Print Survfit Summary
.DN
Prints the result of `summary.survfit'.
.CS
print.summary.survfit(x, digits = max(options() $digits-4, 3), ...)
.RA
.AG x
an object of class `"summary.survfit"', which is the result of the
`summary.survfit' function.
.OA
.AG digits
the number of digits to use in printing the numbers.
.RT
`x', with the invisible flag set to prevent printing.
.SE
prints the summary created by `summary.survfit'.
.SA
`options', `print', `summary.survfit'.
.KW print
.WR
