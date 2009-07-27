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
        # Add parens to the random formula
        rtemp <- formula(paste('(', paste(deparse(random[[2]]), collapse=''), 
                                        ')'))  
        formula[[j]] <- call('+', formula[[j]], rtemp)  # paste it on
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
        if (missing(variance)) theta <- NULL
        else  theta <- variance 
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
            dropx <- c(dropx, temp$terms)
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
        varname <- allnames(x)
        m2 <-  mf[,varname]
        for (i in 1:ncol(m2)) {
            if (is.factor(m2[[i]])) {
                nlev <- length(levels(m2[[i]]))
                contrasts(m2[[i]], nlev) <- diag(nlev)
                }
            }
        model.matrix(terms(x), m2)
        }
    allnames <- function(x){
        if (is.call(x)) {
            if (length(x) == 3) return( c(allnames(x[[2]]), allnames(x[[3]])))
            else return(allnames(x[[2]]))
          }
        if (is.name(x)) as.character(x) else NULL
      }

    getgroups <- function(x, mf) {
        varname <- allnames(x)
        if (is.null(varname)) return(NULL)  # a shrinkage effect like (x1+x2 | 1)
        data.frame(lapply(mf[varname], as.factor))
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
        for (i in 1:nrandom) varlist[[i]] <- coxvarFull #default
        }
    else {
        if (nrandom==1) { # allow a single non-list
            if (!is.list(varlist)) varlist <- list(varlist)
            }
        if (length(varlist) != nrandom) stop ("Wrong length for varlist")
        for (i in 1:length(varlist)) {
            if (!inherits(varlist[[i]], 'coxvar'))
                varlist[[i]] <- coxvarMlist(varlist[[i]])
            }
        }
    fmat <- zmat <- NULL
    nfac <- nslope <- integer(nrandom)
    stemp <- sparse
    for (i in 1:nrandom) {
        f2 <- formula2(flist$random[[i]])
        vfun <- varlist[[i]]
        if (!is.null(f2$interaction)) stop("Interactions not yet written")

        cmat <- getcmat(f2$fixed, m)
        groups <- getgroups(f2$group, m)
        init <- vfun$init(vinit[[i]], variance[[i]], intercept=f2$intercept, 
                            groups, cmat, stemp)
        vparms[[i]] <- init$parms

        if (f2$intercept) {
            if (!is.matrix(init$F) || nrow(init$F) !=n) 
                stop("Invalid result from coxvar function for F")
            nfac[i] <- ncol(init$F)
            fmat <- cbind(fmat, init$F)
            if (stemp[2] < 1) {
                nsparse <- init$sparse
                stemp[2] <- 1
                }
            }
        if (!is.null(cmat)) {
            temp <- newzmat(cmat, fmat)
            zmat <- cbind(zmat, temp)
            nslope[i] <- ncol(temp)
            }
    } 
    if (nsparse>0 & nfac[1]==0) { #must reorder
        # Fix this later
        stop("Only the first random term can be sparse")
        }
    browser()
    fit <- coxme.fit(X, Y, strats, offset, init, control, weights=weights,
                     ties=ties, row.names(m), refine.n,
                     varlist, vparm, thetalist, fixed, 
                     fmat, zmat, refine.n)
    if (is.character(fit)) {
        fit <- list(fail=fit)
        oldClass(fit) <- 'coxme'
        return(fit)
        }
    time2 <- proc.time()
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

    if (ncluster==1) {
        names(fit$frail) <- dimnames(varlist[[1]][[1]])[[1]]
        flinear <- fit$frail[kindex]
        }
    else {
        ftemp <- vector('list', ncluster)
        j <- 0
        flinear <- 0
        for (i in 1:ncluster) {
            tname <- dimnames(varlist[[i]][[1]])[[1]]
            nf <- length(tname)
            temp <- fit$frail[j + 1:nf]
            flinear <- flinear + temp[kindex[,i]]
            names(temp) <- tname
            ftemp[[i]]<- temp
            j <- j+nf
            }
        names(ftemp) <- gnames
        fit$frail <- ftemp
        }

    if (nvar ==0) fit$linear.predictor <- as.vector(flinear)
    else fit$linear.predictor <- as.vector(flinear + c(X %*% fit$coef$fixed))
    fit$n <- nrow(Y)
    fit$terms <- Terms
    fit$assign <- attr(X, 'assign')

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    if (x)  {
        fit$x <- X
        if (length(strats)) fit$strata <- strata.keep
        }
    if (y)     fit$y <- Y
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    time3 <- proc.time()
    timeused <- c((time1[1]+ time1[2]) - (time0[1] + time0[2]), fit$timeused,
                  (time3[1]+ time3[2]) - (time2[1] + time2[2]))
    timeused <- c(sum(timeused), timeused)
    names(timeused) <- c("Total", "setup", "fit1", "fit2", "fit3", "finish")
    fit$timeused <- timeused

    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- Call
    fit$ties <- ties
    fit$kindex <- kindex
    names(fit$loglik) <- c("NULL", "Integrated", "Penalized")
    oldClass(fit) <- 'coxme'
    fit
    }
