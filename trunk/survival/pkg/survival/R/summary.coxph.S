# $Date: 2007-05-10 18:45:52 $ $Id$
summary.coxph <-
 function(object, conf.int = 0.95, scale = 1,
          digits = max(options()$digits - 4, 3)) {
    cox <- object
    if (is.null(cox$coef)) {   # Null model
	return(object)  #The summary method is the same as print in this case
	}

    #
    # Create the matrix of coefficients
    beta <- cox$coef
    nabeta <- !(is.na(beta))          #non-missing coefs
    beta2 <- beta[nabeta]
    if(is.null(beta) | is.null(cox$var))
        stop("Input is not a valid coxph fit")
    se <- sqrt(diag(cox$var))
    if (!is.null(cox$naive.var)) nse <- sqrt(diag(cox$naive.var))

    if (is.null(cox$naive.var)) {
        coef <- cbind(beta, exp(beta), se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
        dimnames(coef) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "z", "p"))
        }
    else {
        coef <- cbind(beta, exp(beta), nse, se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
        dimnames(coef) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "robust se", "z", "p"))
        }

    # Create the matrix of confidence intervals
    z <- qnorm((1 + conf.int)/2, 0, 1)
    beta <- beta * scale
    se <- se * scale
    cimat <- cbind(exp(beta), exp(.Uminus(beta)), exp(beta - z * se),
            exp(beta + z * se))
    dimnames(cimat) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))


    result <- list(coefficient=coef, cimat=cimat)
    result <- c(result, 
                cox[match(c('na.action', 'call', 'fail', 'icc', 'loglik',
                          'score', 'rscore', 'naive.var', 'n'), names(cox),
                        nomatch=0)])
    oldClass(result) <- 'summary.coxph'
    result
    }
