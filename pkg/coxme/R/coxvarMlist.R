# Automatically generated from all.nw using noweb
coxvarMlist <- function(..., rescale=TRUE, pdcheck=TRUE,  positive=NULL) {
    varlist <- list(...)
    # Because of environments, the init function will inherit the
    #  three variables below
    rescale <- rescale
    pdcheck <- pdcheck
    positive <- positive
    
    initialize <- function(initial, fixed, intercept, G, Z, sparse) {
        ncluster <- length(G)
        if (ncluster==0) stop ("Mlist variance requires a grouping variable")
        if (length(Z)>0) stop ("Mlist variance does not allow random slopes")
        if (!intercept)  stop ("Mlist variance applies only to intercepts")

        groups <- expand.nested(G)

        if (length(varlist)==1) varlist <- varlist[[1]] #unlist singletons
        temp <- coxme.varcheck(ncluster, varlist, n=length(G[[1]]),
                               gvars= names(G), 
                               groups= groups[[ncluster]], sparse,
                               rescale, pdcheck)
        ntheta <- temp$ntheta
        theta <- seq(.2, .3, length=ntheta) 

        if (length(initial)>0) {
            if (length(initial) != ntheta) 
                return(list(error="Wrong length for initial vector"))
            theta <- initial
            }
        if (length(fixed) >0) {
            if (length(fixed) != ntheta)
                return(list(error="Wrong length for fixed values"))
            which.fixed <- (!(is.na(fixed) | fixed==0))
            theta[which.fixed] <- fixed
            }
        else which.fixed <- rep(FALSE, ntheta)

        if (is.null(positive)) positive <- rep(TRUE, ntheta)
        else {
            if (!is.logical(positive))
                return(list(
                       error="Positivity constraint must be a logical vector"))
            if (length(positive) != ntheta) 
                return(list(error="Wrong length for positivity constraint"))
            }
        if (any(positive & theta <=0))
            return(list(error="Invalid initial value, must be positive"))        
        theta[positive] <- log(theta[positive])

        list(F=temp$kindex, X=NULL, theta=theta[!which.fixed],
             parms=list(varlist =temp$varlist[[1]], theta=theta,
                        fixed=which.fixed, positive=positive))
        }
    
     generate <- function(newtheta, parms) {
         theta <- parms$theta
         theta[!parms$fixed] <- newtheta
         theta[parms$positive] <- exp(theta[parms$positive])
         
         varmat <- parms$varlist[[1]] * theta[1]
         if (length(theta) >1) {
             for (i in 2:length(theta)) {
                 varmat <- varmat + theta[i]*parms$varlist[[i]]
                 }
             }
         varmat
         }

    wrapup <- function(newtheta,parms) {
        theta <- parms$theta
        theta[!parms$fixed] <- newtheta
        theta[parms$positive] <- exp(theta[parms$positive])
        theta
        }
    
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    oldClass(out) <- 'coxvar'
    out
    }
