# $Date: 2006-08-28 14:16:20 $ $Id$
plot.coxph <- function(fit, ...) {
    op <- par(ask=TRUE)
    on.exit(par(op))
    yy <- (1-fit$residuals)/ fit$linear.predictors   # psuedo y
    plot(fit$linear.predictors, rank(yy))

    std.resid <- fit$residuals/ sqrt(predict(fit, type='expected'))
    plot(fit$linear.predictors, std.resid,
	xlab='Linear predictor', ylab='standardized residual')

    }
