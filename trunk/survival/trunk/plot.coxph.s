#SCCS $Date: 1992-04-14 18:07:07 $ $Id: plot.coxph.s,v 4.2 1992-04-14 18:07:07 grill Exp $
plot.coxph <- function(fit, ...) {
    op <- par(ask=T)
    on.exit(par(op))
    yy <- (1-fit$residuals)/ fit$linear.predictors   # psuedo y
    plot(fit$linear.predictors, rank(yy))

    std.resid <- fit$residuals/ sqrt(predict(fit, type='expected'))
    plot(fit$linear.predictors, std.resid,
	xlab='Linear predictor', ylab='standardized residual')

    }
