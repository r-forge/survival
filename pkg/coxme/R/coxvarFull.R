# Automatically generated from all.nw using noweb
coxvarFull <- function(collapse=FALSE) {
    collapse <- collapse
    # Because of R's lexical scoping, the values of the options
    #  above, at the time the line below is run, are known to the
    #  initialize function
    initialize <- function(initial, fixed, intercept, G, Z,  sparse) {
        if (is.null(Z)) nvar <- 0
        else nvar <- ncol(Z)
        initmatch <- function(namelist, init) {
            if (is.null(names(init))) iname <- rep('', length(init))
            else iname <- names(init)
            
            indx1 <- pmatch(iname, namelist, nomatch=0, duplicates.ok=TRUE)
            if (any(iname=='')) {
                temp <- 1:length(namelist)
                if (any(indx1>0)) temp <- temp[-indx1]   #already used
                indx2 <- which(iname=='')
                n <- min(length(indx2), length(temp))
                if (n>0) indx1[indx2[1:n]] <- temp[1:n]
                }
            indx1
            }
        
        if (intercept & is.null(G))
            return(list(error=("Invalid null random term (1|1)")))
        
        if (intercept && nvar==0) {
            gname <- names(G)
            ntheta <- length(gname)
            theta <- rep(.2, ntheta)
            if (length(initial) >0) {
                temp <- initmatch(gname, initial)
                if (any(temp==0))
                    return(list(error=paste('Element', which(temp==0),
                                            'of initial values not matched')))
                else theta[temp] <- unlist(initial)
                if (any(theta <=0))
                    return(list(error='Invalid initial value'))
                }

            which.fixed <- rep(FALSE, ntheta)
            if (length(fixed)>0) {
                temp <- initmatch(gname, fixed)
                if (any(temp==0))
                    return(list(error=paste('Element', which(temp==0),
                                             'of variance values not matched')))
                else theta[temp] <- unlist(fixed)
                which.fixed[temp] <- TRUE
                }
            if (ncol(G) ==1) {
                gtemp <- as.factor(G[[1]])[,drop=TRUE] #drop unused levels
                nlevel <- length(levels(gtemp))
                gfrac <- table(gtemp)/ length(gtemp)
                if (nlevel >= sparse[1] && any(gfrac <= sparse[2])) {
                    indx <- order((gfrac> sparse[2]), 1:nlevel)  #False then True for order
                    gtemp <- factor(gtemp, levels=levels(gtemp)[indx])
                    nsparse <- sum(gfrac <= sparse[2])
                    if (nsparse== nlevel) vmat<- bdsI(nsparse)
                    else {
                        k <- nlevel - nsparse  #number of non-sparse levels
                        rmat <- matrix(0., nrow=nlevel, ncol=k)
                        rmat[seq(nsparse+1, by= nlevel+1, length=k] <- 1.0
                        vmat <- bdsmatrix(blocksize=rep(1,nsparse), 
                                          blocks= rep(1,nsparse), rmat=rmat)
                        }
                    }
                else vmat <- diag(nlevel)
                list(F=matrix(as.numeric(gtemp)), X=NULL, 
                     theta=log(theta[!which.fixed]), 
                     parms=list(vmat=vmat, theta=theta, level=levels(gtemp),
                                fixed=which.fixed, case=1, gname=gname))
                    }
            else {
                if (!collapse) {
                    G <- expand.nested(G)
                    n.nest <- ncol(G)
                    F <- matrix(0, nrow=nrow(G), ncol=n.nest)
                    nlevel <- sapply(G, function(x) length(levels(x)))
                    levellist <- vector('list', n.nest)
                    
                    # Sparsity?
                    gtemp <- G[,n.nest]
                    gfrac <- table(gtemp)/ length(gtemp)
                    if (nlevel[ncol(G)] > sparse[1] && any(gfrac <= sparse[2])) {
                        indx <- order((gfrac> sparse[2]), 1:nlevel)
                        gtemp <- factor(gtemp, levels=levels(gtemp)[indx])
                        nsparse <- sum(gfrac <= sparse[2])
                    
                        F[,1] <- as.integer(gtemp) 
                        levellist[[1]] <- levels(gtemp)
                        nlevel <- rev(nlevel)
                        gname <- rev(gname)
                        theta <- rev(theta)
                        which.fixed <- rev(which.fixed)
                        
                        for (i in 2:n.nest) {
                            j <- 1 + n.nest -i
                            F[,i] <- as.numeric(G[,j])
                            levellist[[i]] <- levels(G[,j])
                            }
                        }
                    else { # No sparse, so don't reverse the levels & annoy the user
                        nsparse <- 0
                        for (i in 1:n.nest) {
                            F[,i] <- as.numeric(G[,i])
                            levellist[[i]] <- levels(G[,i])
                            }
                        }

                    list(F=F, X=NULL, 
                             theta=log(theta[!which.fixed]),
                             parms=list(nlevel=nlevel, nsparse=nsparse, gname=gname,
                                    fixed=which.fixed, levels=levellist, 
                                    theta=theta, case=2))
                    }
                else { #ncol(G)>1, intercept=T, nvar=0, collapse=T
                    gtemp <- as.factor(G[,1])[,drop=TRUE]
                    gfrac <- table(gtemp)/ n
                    if (nlevel[ncol(G)] > sparse[1] && any(gfrac <= sparse[2])) {
                        indx <- order((gfrac> sparse[2]), 1:nlevel)
                        G[,1] <- factor(gtemp, levels=levels(gtemp)[indx])
                        nsparse <- sum(gfrac <= sparse[2])
                        }
                    else nsparse <- 0  
                    
                    G <- expand.nested(G)
                    indx <- unique(G[,ncol(G)])
                    ncoef <- length(indx)
                    varlist <- vector('list', ncol(G))
                    if (ncol(G)>2) {for (i in 2:(ncol(G)-1)) 
                        varlist[[i]] <- bdsBlock(1:ncoef, G[indx,i])
                          }
                    varlist[ncol(G)] <- bdsI(ncoef)
                    
                    temp <- bdsBlock(1:ncoef, G[indx,1])
                    if (nsparse>0 && all(gfrac <= sparse[2])) varlist[[1]] <- temp
                    else { #pick off part
                        temp <- bdsBlock(1:ncoef, G[indx,1])
                        tsize <- temp@blocksize[1:max(1,nsparse)]
                        dense <- (1 + sum(tsize)):ncoef
                        rtemp <- matrix(0, sum(tsize), ncoef)
                        rmat[by=nrow(rmat)+1, to=length(rmat), length=ncol(rmat)] <- 1.0
                        varlist[[1]] <- bdsmatrix(blocksize=tsize,
                                              blocks=temp@blocks[1:sum(tsize*(tsize+1)/2)],
                                              rmat=rbind(rtemp, as.matrix(temp[dense,dense])))
                        }
                    
                    for (i in 2:ncol(G)) varlist[[i]] <- varlist[[i]] + 0*varlist[[1]]
                      
                    list(F=matrix(as.numeric(ulist)), X=NULL, 
                         theta=log(theta[!which.fixed]),
                         parms=list(varlist=varlist, theta=theta, 
                                    fixed=which.fixed, gname=gname,
                                    levels=lapply(G, function(x) levels(x)), case=2.5))
                      }
                    }
            }
        else if (is.null(G)) {
            xname <- dimnames(Z)[[2]]
            if (length(initial) >0) {
              temp <- initmatch(xname[1], initial)
              if (any(temp==0)) 
                  return(list(error=paste('Element', which(temp==0),
                                          'of initial values not matched')))
              else theta <- initial
              }
            else theta <- .2 / mean(sqrt(apply(Z,2,var)))
              
            if (length(fixed) >0) {
                temp <- initmatch(xname[1], fixed)
                if (any(temp==0))
                    return(list(error=paste('Element', which(temp==0),
                                            'of fixed variance values not matched')))
                else theta <- fixed
                which.fixed <- TRUE
                }
            else which.fixed <- FALSE
            if (theta <=0) return(list(error="Invalid variance value, must be >0"))

            list(theta=theta[!which.fixed], F=NULL, X=Z,
                     parms=list(fixed=which.fixed, theta=theta,
                                xname=xname, case=3))
                }
                }
        else {
            gtemp <- as.factor(G[[1]])[,drop=TRUE] #drop unused levels
            nlevel <- length(levels(gtemp))
            gfrac <- table(gtemp)/ length(gtemp)
            if (nlevel[1] > sparse[1] && any(gfrac <= sparse[2])) {
                indx <- order((gfrac> sparse[2]), 1:nlevel)
                G[,1] <- factor(gtemp, levels=levels(gtemp)[indx])
                nsparse <- sum(gfrac <= sparse[2])
               }
            else nsparse <- 0

            G <- expand.nested(G)
            ngroup <- ncol(G)
            F <- matrix(0, nrow=nrow(G), ncol=n.nest)
            nlevel <- sapply(G, function(x) length(levels(x)))
            levellist <- lapply(G, levels)
            for (i in 1:n.nest) 
                F[,i] <- as.numeric(G[,i])

            if (ngroup==1 && nsparse==nlevel[1]) vmat <- bdsI(nlevel[1])
            else {
                if (nsparse==0) { #make the whole first term a block
                    rmat <- matrix(0., nrow=sum(nlevel), ncol=sum(nlevel[-1]))
                    rmat[nlevel[1]+1, by=nrow(rmat) +1, length=ncol(rmat)] <- 1.0
                    vmat <- bdsmatrix(blocksize=nlevel[1], blocks=diag(nlevel[1]),rmat=rmat)
                    }
                else {
                    rmat <- matrix(0., nrow=sum(nlevel), ncol=sum(nlevel)-nsparse)        
                    rmat[nlevel[1]+1, by=nrow(rmat) +1, length=ncol(rmat)] <- 1.0
                    vmat <- bdsmatrix(blocksize=rep(1,nsparse), blocks=rep(1., nsparse),
                                      rmat=rmat)
                    }
            X <- matrix(0., nrow=n, ncol=sum(nlevel)*nvar)
            indx <- seq(0, length=nvar, by=sum(nlevel))
            offset <-0
            for (i in 1:ngroup) { 
                for (j in 1:nvar) {
                    for (k in 1:nvlevel[i])
                        X[,offset+ k +indx] <- Z[,j] * (F[,i]==k)
                    offset <- offset + nlevel[i]
                    }
                }
            }
        }
       generate= function(newtheta, parms) {
           theta <- parms$theta
           if (length(newtheta)>0) theta[!parms$fixed] <- exp(newtheta)

           if (parms$case==1) return(theta*parms$vmat)
    if (parms$case==2) {
        temp <- rep(theta, parms$nlevel)
        if (parms$nsparse >0) {
            bdsmatrix(blocksize=rep(1, parms$nsparse),
                      blocks=temp[1:parms$nsparse],
                      rmat=rbind(matrix(0.,nrow=parms$nsparse, ncol=sum(nlevel)),
                                 diag(temp[-(1:parms$nsparse)])))
            }
         else diag(temp)
         }
    if (parms$case== 2.5) {
        temp <- parms$varlist[[1]] * theta[1]
        for (i in 2:length(theta))
            temp <- temp + parms$varlist[[i]] * theta[i]
        return(temp)
        }
    if (parms$case==3) return(diag(length(parms$xname)) * theta)
    if (parms$case==4) {
        ngroup <- parms$ngroup
        nvar <- parms$nvar

        vmat <- diag(nvar+ngroup)  #var/cov matrix = transformed parameters
        temp <- exp(parms)            
        std <- sqrt(temp[1:(nvar+ngroup)])
        offset <- ngroup + nvar

        for (i in 1:ngroup) {
            temp2 <- temp[1:nvar +offset]
            vmat[i, -(1:ngroup)] <- (temp2-1)/(temp2+1)
            }
        for (i in 1:nvar) {
            if (i<nvar){
                temp2 <- temp[1:(nvar-i) + offset]
                vmat[i+ngroup, i+ngroup +1:(nvar-i)] <- (temp2-1)/(temp2+1)
                offset <- offset + nvar -i
                }
            }   
        vmat <- diag(std) %*% vmat %*% diag(std)   
        n1 <- sum(parms$nlevel)  #number of intercept coefs
        n2 <- n1*nvar            #number of slope coefs
        rcol <- ncol(parms$vmat@rmat)
        rtemp <- matrix(0., nrow=n1 + n2, ncol=n2+ rcol)
        indx1 <- seq(1, by=nrow(rtemp)+1, length=ngroup)
        indx2 <- 0   
        for (i in 1:ngroup) {
            for (j in 1:nvar) {
                indx3 <- (rcol+j-1)*nrow(rtemp)
                rtemp[indx1+indx2+indx3] <- vmat[i,j]
                }
            indx2 <- indx2 + parms$nlevel[i]
            } 
        if (nvar >1) {
            for (i in 1:ngroup) {
                for (j in 1:(nvar-1)) {
                    for (k in (j+1):nvar) {
                        indx3 <- (rcol+j-1)*nrow(rtemp)
                        rtemp[indx1+indx2+indx3] <- vmat[ngroup+i, ngroup+j]
                        }
                    }
                indx2 <- indx2 + parms$nlevel[i]
                }
            }

        final <- bdsmatrix(blocksize=parms$vmat@blocksize,
                          blocks=parms@vat@blocks,
                          rmat=rtemp)
        diag(final) <- rep(diag(vmat), rep(parms$nlevel, nvar+1))
        final
        }
    wrapup <- function(theta, b, parms) {
        if (parms$case <4) {
            newtheta <- parms$theta
            newtheta[!parms$fixed] <- exp(theta)
            }
        
        if (parms$case==1) {
            names(theta) <- parms$gname
            names(b) <- parms$levels
            return(list(theta=newtheta, random=b))
            }

        if (parms$case==2) {
            ngroup <- length(parms$nlevel)
            names(b) <- unlist(parms$levellist)
            random <- split(b, rep(1:ngroup, parms$nlevel))
            names(random) <- parms$tname
            return(list(theta=newtheta, random=random))
            }

        if (parms$case==3) {
            names(theta) <- parms$xname
            names(b) <- parms$levels
            return(list(theta=newtheta, random=b))
            }
        
        if (parms$case==4) {
             ngroup <- length(parms$nlevel)
             nvar <- parms$nvar
             names(b) <- rep(unlist(parms$levellist, 1+nvar))
             random <- split(b, rep(rep(1:ngroup, parms$nlevel), 1+nvar))
             names(random) <- parms$tname[1:(ngroup+nvar)]
             
             cmat <- diag(nvar+ngroup)  #correlation matrix
             temp <- exp(parms$theta)
             temp[!parms$fixed] <- exp(theta)
             offset <- ngroup + nvar

             for (i in 1:ngroup) {
                 temp2 <- temp[1:nvar +offset]
                 cmat[i, -(1:ngroup)] <- (temp2-1)/(temp2+1)
                 }
             for (i in 1:nvar) {
                 if (i<nvar){
                     temp2 <- temp[1:(nvar-i) + offset]
                     cmat[i+ngroup, i+ngroup +1:(nvar-i)] <- (temp2-1)/(temp2+1)
                     offset <- offset + nvar -i
                     }
                 }   
             cmat <- cmat + t(cmat)
             diag(cmat) <- 1.0
             
             indx <- 1:(nvar+ngroup)
             vars <- temp[]
             names(vars) <- temp[indx]
             names(vars) <- parms$tname[indx]
             dimnames(cmat) <- list(parms$tname, parms$tname)
             return(list(theta=list(variance=vars, correlation=cmat),
                  random=random))
             }
        
        if (parms$case==2.5) stop("Wrapup not finished")
        }
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    oldClass(out) <- 'coxvar'
    out
    }
