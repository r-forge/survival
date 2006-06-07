#
# Create a data set similar to the one provided by Al Zinsmeister
#  It is a hard test case for survfit.turnbull
#
time1 <- c(rep(0,100), rep(1,200), 100, 200, 210, 220,
           rep(365,100), rep(366,5), 731:741)

time2 <- c((1:100)*3,  10+1:100, rep(365:366, c(60,40)), NA, 500, NA, 450,
           rep(730,90), rep(NA,10), c(528,571,691,730,731),
           NA, 1095:1099, NA, 1400, 1200, 772, 1461)

zfit <- survfit(Surv(time1, time2, type='interval2') ~1)

#
# There are 100 intervals of the form (0,x) where x is from 3 to 300,
#  and 200 more of the form (1,x) where x is from 11 to 366.  These
#  lead to a mass point in the interval (1,3), which is placed at 2.
#  The starting estimate has far too little mass placed here, and it takes
#  the EM a long time to realize that most of the weight for the first 300
#  subjects goes here.  With acceleration, it takes 16 iterations, without
#  it takes >40.  (On Al's orginal data, without accel still wasn't there after
#  165 iters!)
#
# The next 4 obs give rise to potential jumps at 100.5, 200.5, 211.5, and
#  221.  However, the final estimate has no mass at all on any of these.
#  Assume mass of a,b, and c at 2, 100.5 and 365.5, and consider the 
#  contributions: 
#    123 obs that overlap a only
#    137 obs that overlap a and b
#     40 obs that overlap a, b, c
#      1 obs that overlap b, c
#    108 obs that overlap c   (200, 210,200, 365, and 366 starting points)
#  For some trial values of a,b,c, compare the loglik to that of (a+b),0,c
#   First one: a^123 (a+b)^137 (a+b+c)^40 (b+c) c^108
#   Second:    (a+b)^123 (a+b)^137 (a+b+c)^40 c c^108
#   Likelhood improves if (1 + b/a)^123 > 1+ b/c, which is true for almost
#     all a and c.  In particular, at the solution a and c are approx .7 and
#     .18, respectively.
#
# The program can't see this coming, of course, and so iterates towards a
#  KM with epsilon sized jumps at 100.5, 200.5, and 211.5.  Whether these
#  intervals should be removed during iteration, as detected, is an open
#  question for me.
#
#
# True solution: mass points at 2, 365.5, 408, and 756.5, of sizes a, b, c, d
# Likelihood:    a^260 (a+b)^40 (b+c)^92 (b+c+d)^12 c^5 d^11
# Solution: a=0.6958, b=0.1674, c=0.1079, d=0.0289

tfun <- function(x) {
    if (length(x) ==3) x <- c(x, .03)
    x <- x/sum(x)  #make probabilities sum to 1
    loglik <- 260*log(x[1]) + 40*log(x[1]+x[2]) + 92*log(x[2] + x[3]) +
                       12*log(x[2]+x[3]+x[4]) + 5*log(x[3]) + 11*log(x[4])
    -loglik  #find the max, not the min
    }

nfit <- nlminb(start=c(.7,.15, .1), tfun, lower=0, upper=1)
nparm <- c(nfit$para, .03)
nparm <- nparm / sum(nparm)
zparm <- -diff(c(1, zfit$surv[match(c(2, 365.5, 408, 756.5), zfit$time)]))
aeq(round(tfun(nparm),4), round(tfun(zparm),4))
# .0001 is the tolerance in survfit.turnbull

rm(tfun, nfit, nparm, zparm, time1, time2, zfit) 
