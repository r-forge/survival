.BG
.FN rowsum
.TL
Give row sums of a matrix, based on a grouping variable.
.CS
rowsum(x, group, reorder=T)
.RA
.AG x
 a matrix or vector of numeric data.  Missing values are allowed.
.AG group
 a vector giving the grouping, with one element per row of `x'.
Missing values are not allowed.
.OA
.AG reorder
if True, then the result will be in order of sort(unique(group)),
if False, it will be in the order that rows were encountered (and
may run faster for large matrices).
The default is to reorder the data, so as to agree with tapply (see
example below).
.RT
a matrix containing the sums.  There will be one row per unique value
of `group'.
.SA
tapply
.EX
x <- matrix(runif(100), ncol=5)
group <- sample(1:8, 20, T)
xsum <- rowsum(x, group)

#same result another way, slower, and temp may be much larger than x
temp <- model.matrix( ~a -1, data.frame(a=as.factor(group)))
xsum2<- t(temp) %*% x

#same as last one, but really slow
xsum3 <- tapply(x, list(group[row(x)], col(x)), sum)

.KW manip
.WR
