#SCCS $Date: 1992-03-04 16:48:11 $ $Id: plot.coxph.s,v 4.1 1992-03-04 16:48:11 therneau Exp $
plot.coxreg <- function(fit, ...) {
    op <- par(ask=T)
    on.exit(par(op))
    yy <- (1-fit$residuals)/ fit$linear.predictors   # psuedo y
    plot(fit$linear.predictors, rank(yy))

    std.resid <- fit$residuals/ sqrt(predict(fit, type='expected'))
    plot(fit$linear.predictors, std.resid,
	xlab='Linear predictor', ylab='standardized residual')

    }
