.BG
.FN rowsum
.TL
Give row sums of a matrix, based on a grouping variable.
.CS
rowsum(x, group)
.RA
.AG x
 a matrix or vector of numeric data.
.AG group
 a vector giving the grouping, with one element per row of `x'.
.RT
a matrix containing the sums.  There will be one row per unique value
of `group'.
.SA
tapply
.EX
x <- matrix(runif(100), ncol=5)
group <- sample(1:8, 20, T)
xsum <- rowsum(x, group)

#same result another way, slower, rows are ordered by sort(unique(group))
temp <- model.matrix( ~a -1, data.frame(a=as.factor(group)))
xsum3<- t(temp) %*% x

#same as last one, but really slow
xsum2 <- tapply(x, list(group[row(x)], col(x)), sum)

.KW manip
.WR
