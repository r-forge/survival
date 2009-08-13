# Automatically generated from all.nw using noweb
coxme <- function(formula,  data, 
        weights, subset, na.action, init, 
        control, ties= c("efron", "breslow", "exact"),
        singular.ok =T, varlist=NULL, variance, vinit=.2, sparse=c(50,.02),
        rescale=TRUE, pdcheck=TRUE, x=FALSE, y=TRUE, shortlabel=TRUE, 
        refine.n=0, random, fixed, ...) {

    time0 <- proc.time()    #debugging line
    ties <- match.arg(ties)
    Call <- match.call()

    if (!missing(fixed)) {
        if (missing(formula)) {
            formula <- fixed
            warning("The 'fixed' argument of coxme is depreciated")
            }
        else stop("Both a fixed and a formula argument are present")
        }
    if (!missing(random)) {
        warning("The random argument of coxme is depreciated")
        if (class(random) != 'formula' || length(random) !=2) 
            stop("Invalid random formula")
        j <- length(formula)   #will be 2 or 3, depending on if there is a y

        # Add parens to the random formula and paste it on
        formula[[j]] <- call('+', formula[[j]], call('(', random[[2]]))  
        }
    temp <- call('model.frame', formula= subbar(formula))
    for (i in c('data', 'subset', 'weights', 'na.action'))
        if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
        Y <- model.extract(m, "response")
        n <- nrow(Y)
        if (!inherits(Y, "Surv")) stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type!='right' && type!='counting')
            stop(paste("Cox model doesn't support '", type,
                              "' survival data", sep=''))

        weights <- model.weights(m)
        if (length(weights) ==0) weights <- rep(1.0, n)
        else if (any(weights <=0))
            stop("Negative or zero weights are not allowed")

        offset <- model.offset(m)
        if (length(offset)==0) offset <- rep(0., n)

        # Check for penalized terms; the most likely is pspline
        pterms <- sapply(m, inherits, 'coxph.penalty')
        if (any(pterms)) {
            stop("You cannot have penalized terms in coxme")
            }

        if (missing(control)) control <- coxme.control(...)
        if (missing(init)) init <- NULL
        flist <- formula1(formula)
        if (hasAbar(flist$fixed))
            stop("Invalid formula: a '|' outside of a valid random effects term")

        special <- c("strata", "cluster")
        Terms <- terms(flist$fixed, special)
        attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
        strats <- attr(Terms, "specials")$strata
        cluster<- attr(Terms, "specials")$cluster
        if (length(cluster)) {
            stop ("A cluster() statement is invalid in coxme")
            }
        if (length(strats)) {
            temp <- untangle.specials(Terms, 'strata', 1)
            if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
            else strata.keep <- strata(m[,temp$vars], shortlabel=T)
            strats <- as.numeric(strata.keep)
            X <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
            }
        else X <- model.matrix(Terms, m)[,-1,drop=F]


    nrandom <- length(flist$random)
    if (nrandom ==0) stop("No random effects terms found")
    getcmat <- function(x, mf) {
        if (is.null(x)) return(NULL)
        Terms <- terms(eval(call("~", x)))
        attr(Terms, 'intercept') <- 0  #ignore any "1+" that is present

        varnames <-  attr(Terms, 'term.labels')
        ftemp <- sapply(mf[varnames], is.factor)
        if (any(ftemp)) {
            clist <- lapply(mf[ftemp], function(x) diag(length(levels(x))))
            model.matrix(Terms, mf, contrasts.arg =clist)
            }
        else model.matrix(Terms, mf)
        }
    getGroupNames <- function(x) {
        if (is.call(x) && x[[1]]==as.name('/')) 
            c(getGroupNames(x[[2]]), getGroupNames(x[[3]]))
        else deparse(x)
        }

    getgroups <- function(x, mf) {
        varname <- getGroupNames(x)
        if (varname=='1') return(NULL)  # a shrinkage effect like (x1+x2 | 1)
        else indx <- match(varname, names(mf), nomatch=0)
        if (any(indx==0)) stop(paste("Invalid grouping factor", varname[indx==0]))
        else data.frame(lapply(mf[indx], as.factor))
        }
    if (missing(vinit)) vinit <- vector('list', nrandom)
    else {
        if (nrandom==1 && is.numeric(vinit)) vinit <- list(vinit)
        if (!is.list(vinit)) stop("Invalid value for `vinit` parameter")
        if (length(vinit) > nrandom) stop ("Invalid length for vinit")
        if (!all(sapply(vinit, is.numeric))) 
            stop("Vinit must contain numeric values") 
        
        if (length(vinit) < nrandom) 
            vinit <- c(vinit, vector('list', nrandom - length(vinit)))
                       
        tname <- names(vinit)
        if (!is.null(tname)) {
            temp <- pmatch(tname, names(flist$random), nomatch=0)
            temp <- c(temp, (1:nrandom)[-temp])
            vinit <- vinit[temp]
            }
      }

    if (missing(variance)) variance <- vector('list', nrandom)
    else {
        if (nrandom==1 && is.numeric(variance)) variance <- list(variance)
        if (!is.list(variance)) stop("Invalid value for `variance` parameter")
        if (length(variance) > nrandom) stop ("Invalid length for variance")
        if (!all(sapply(variance, is.numeric))) 
            stop("Variance must contain numeric values") 
        
        if (length(variance) < nrandom) 
            variance <- c(variance, vector('list', nrandom - length(variance)))
                       
        tname <- names(variance)
        if (!is.null(tname)) {
            temp <- pmatch(tname, names(flist$random), nomatch=0)
            temp <- c(temp, (1:nrandom)[-temp])
            variance <- variance[temp]
            }
      }
    newzmat <- function(cmat, fmat) {
        if (length(fmat)==0) return(cmat)  # a formula with (... | 1)
        newcol <- ncol(cmat) * sum(apply(fmat,2,max))
        newz <- matrix(0., nrow=nrow(cmat), ncol=newcol)
        indx <- 0
        nc <- ncol(cmat)
        for (i in 1:ncol(fmat)){
            for (j in 1:max(fmat[,i])) {
                newz[, 1:nc + indx] <- cmat * (fmat[,i]==j)
                indx <- indx + nc
                }
            }
        newz
        }
    vparms <- vector('list', nrandom)
    fname <- zname <- thetalist <- vparms
    if (missing(varlist)) {
        varlist <- vector('list', nrandom)
        for (i in 1:nrandom) varlist[[i]] <- coxvarFull() #default
        }
    else {
        if (nrandom==1) { # allow a single non-list
            if (inherits(varlist, 'coxvar') || !is.list(varlist)) 
                varlist <- list(varlist)
            }
        if (length(varlist) != nrandom) stop ("Wrong length for varlist")
        for (i in 1:length(varlist)) {
            if (!inherits(varlist[[i]], 'coxvar'))
                varlist[[i]] <- coxvarMlist(varlist[[i]])
            }
        }
    fmat <- zmat <- matrix(0, nrow=n, ncol=0)
    ntheta <- integer(nrandom)
    theta <-  NULL   #initial values of parameters to iterate over
    for (i in 1:nrandom) {
        f2 <- formula2(flist$random[[i]])
        if (f2$intercept & f2$group==1)
            stop(paste("Error in random term ", i, 
                       ": Random intercepts require a grouping variable", sep=''))
        vfun <- varlist[[i]]
        if (!is.null(f2$interaction)) stop("Interactions not yet written")

        cmat <- getcmat(f2$fixed, m)
        groups <- getgroups(f2$group, m)
        ifun <- vfun$initialize(vinit[[i]], variance[[i]], intercept=f2$intercept, 
                            groups, cmat, sparse)
        if (!is.null(ifun$error)) 
            stop(paste("In random term ", i, ": ", ifun$error, sep=''))
        vparms[[i]] <- ifun$parms

        theta <- c(theta, ifun$theta)
        ntheta[i] <- length(ifun$theta)

        if (f2$intercept) {
            if (!is.matrix(ifun$F) || nrow(ifun$F) !=n) 
                stop(paste("In random term ", i, 
                           ": Invalid intercept matrix F", sep=''))
            for (i in 1:ncol(ifun$F)) {
                temp <- as.integer(factor(ifun$F[,i]))
                if (any(temp != ifun$F[,i]))
                    stop(paste("In random term ", i,
                               ": intercept matrix has an invalid column", sep=''))
                }
            fmat <- cbind(fmat, ifun$F)
            }

        if (!is.null(cmat)) {
            temp <- newzmat(cmat, fmat)
            zmat <- cbind(zmat, temp)
            }
        } 
    fit <- coxme.fit(X, Y, strats, offset, init, control, weights=weights,
                     ties=ties, row.names(m),
                     fmat, zmat, varlist, vparms, 
                     theta, ntheta, refine.n)
    if (is.character(fit)) {
        fit <- list(fail=fit)
        oldClass(fit) <- 'coxme'
        return(fit)
        }
    fcoef <- fit$coefficients$fixed
    nvar <- length(fcoef)
    if (length(fcoef)>0 && any(is.na(fcoef))) {
        vars <- (1:length(fcoef))[is.na(fcoef)]
        msg <-paste("X matrix deemed to be singular; variable",
                        paste(vars, collapse=" "))
        if (singular.ok) warning(msg)
        else             stop(msg)
        }
    if (length(fcoef) >0) {
        names(fcoef) <- dimnames(X)[[2]]
        fit$coefficients <- list(fixed=fcoef, random=fit$coeff$random)
        }
    rlinear <- rep(0., nrow(Y))
    indx <- 0
    if (length(fmat)) {
        for (i in 1:ncol(fmat)) {
            rlinear <- rlinear + fit$frail[fmat[,i]+indx]
            indx <- indx + max(fmat[,i])
            }
        }
    if (length(zmat)) {
        for (i in 1:ncol(zmat))
            rlinear <- rlinear + fit$frail[indx+i]*zmat[,i]
        }

    if (nvar==0) fit$linear.predictor <- rlinear
    else fit$linear.predictor <- as.vector(rlinear + c(X %*% fit$coef$fixed))
    fit$n <- nrow(Y)
    fit$terms <- Terms
    fit$assign <- attr(X, 'assign')
    fit$formulaList <- flist

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    if (x)  {
        fit$x <- X
        if (length(strats)) fit$strata <- strata.keep
        }
    if (y)     fit$y <- Y
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- Call
    fit$ties <- ties
    names(fit$loglik) <- c("NULL", "Integrated", "Penalized")
    oldClass(fit) <- 'coxme'
    fit
    }
