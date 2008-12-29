survfit <- function (formula, ...) {
    Call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it. 
    # At one time, we allowed this to act like a formula, by adding the ~1
    #    on for the user.  Very non-standard, and a bad idea.
    # We removed all refernce to this useage from the documentation several
    #    years ago.  The line below is a temporary, "nice" failure message
    #    for this no-longer supported use.  In a couple more years, axe this
    #    section.
    #
    if ((mode(Call[[2]]) == 'call' &&  Call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
	stop(paste("Survfit requires a formula or a coxph fit as the",
		   "first argument"))
	}

    UseMethod('survfit', formula)
    }
