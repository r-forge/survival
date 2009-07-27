# Automatically generated from all.nw using noweb

 coxme.fit <- function(x, y, strata, offset, control,
                         weights, ties, rownames, 
                         fmat, zmat, varlist, theta, ntheta,
                         refine.n) {
     time0 <- proc.time()

     n <-  nrow(y)
     if (length(x) ==0) nvar <-0
     else nvar <- ncol(as.matrix(x))
     
     if (missing(offset) || is.null(offset)) offset <- rep(0.0,n)
     if (missing(weights)|| is.null(weights))weights<- rep(1.0,n)
     else {
         if (any(weights<=0)) stop("Invalid weights, must be >0")
         }
     if (ncol(y) ==3) {
         if (length(strata) ==0) {
             sorted <- cbind(order(-y[,2], y[,3]), 
                             order(-y[,1]))
             newstrat <- n
             }
         else {
             sorted <- cbind(order(strata, -y[,2], y[,3]),
                             order(strata, -y[,1]))
             newstrat  <- cumsum(table(strata))
             }
         status <- y[,3]
         ofile <-  'agfit6b'
         rfile <-  'agfit6d'
         }
     else {
         if (length(strata) ==0) {
             sorted <- order(-y[,1], y[,2])
             newstrat <- n
             }
         else {
             sorted <- order(strata, -y[,1], y[,2])
             strata <- (as.numeric(strata))[sorted]
             newstrat <-  cumsum(table(strata))
             }
         status <- y[,2]
         ofile <- 'coxfit6b' # fitting routine
         rfile <- 'coxfit6d' # refine.n routine
         }
     if (!is.list(varlist)) stop("variance matrix list isn't a list!")
     if (is.matrix(fmat) && ncol(fmat >1)) {
         ncluster <- ncol(fmat)
         clnames <- dimnames(fmat)[[2]]
         }
     else ncluster <- 1
     if (ncluster != length(varlist))
         stop("Lengths of variance list and of fmat disagree")

     #
     # Check out fmat:
     #   each column should be integer
     #   each col must contain numbers from 1 to k for some k
     #   column 1 should contain at least one "1" (the C routine
     if (any(fmat != floor(fmat))) stop("fmat must be integers")
     if (any(fmat <1)) stop("fmat must be >0")
     fmat <- as.matrix(fmat)
     nfrail <- apply(fmat, 2, max)
     temp <- apply(as.matrix(fmat), 2, function(x) length(unique(x)))
     if (any(nfrail != temp)) 
         stop("fmat must be a set of integer indices, with none missing")    
  dummy <- kfun(theta, varlist)
  if (is.null(dummy@rmat)) rtemp <- 0
      else                 rtemp <- ncol(dummy@rmat)

  for (i in 2:ncol(fmat)) fmat[,i] <- fmat[,i] + max(fmat[,i-1])
  findex <- matrix(0, nrow=max(fmat), ncol=ncol(fmat))
  for (i in 1:ncol(fmat)) findex[cbind(fmat[,i], i)] <- 1
      
  ifit <- .C('coxfit6a', 
                 as.integer(n),
                 as.integer(nvar),
                 as.integer(ncol(y)),
                 as.double(c(y)),
                 as.double(x),
                 as.double(offset),
                 as.double(weights),
                 as.integer(length(newstrat)),
                 as.integer(newstrat),
                 as.integer(sorted-1),
                 as.integer(ncol(fmat)),
                 as.integer(fmat-1),
                 as.integer(findex),
                 as.integer(length(dummy@blocksize)),
                 as.integer(dummy@blocksize),
                 as.integer(rtemp),
                 means = double(temp.nvar),
                 scale = double(temp.nvar),
                 as.integer(ties=='efron'),
                 as.double(control$toler.chol),
                 as.double(control$eps),
                 as.integer(control$sparse.calc))
      means   <- ifit$means
      scale   <- ifit$scale
 fit0 <- coxph(y ~ x + offset(offset), weights=weights, method=ties)
 whichterm <- rep(1:nrandom, nfac)  # which term does each column of fmat go to
 kfun <- function(theta, varlist, parmlist, 
                  fcount=tapply(nfac.nevel,whichterm, sum), 
                  zcount = nslope) {

     nrandom <- length(varlist)
     if (nrandom == 1) return(varlist[[1]]$generate(theta, parmlist[[1]]))
     # Need to build up the matrix by pasting up a composite R
     nrow.R <- sum(fcount) + sum(zcount)
     ncol.R <- nrow.R - nsparse
     R <- matrix(0., nrow.R, ncol.R)
     
     fcount2 <- fcount; fcount2[1] <- fcount[1] - nsparse  
     indx1 <- cumsum(c(0, fcount2))  #offsets for intercept columns
     indx2 <- cumsum(c(0, fcount))   #offsets for intercept  rows
     indx3 <- cumsum(c(sum(fcount2), zcount)) # offsets for slope cols
     indx4 <- indx3 + nsparse  #     #offsets for slope rows
     for (i in 1:nrandom) {
         temp <- varlist[[i]]$generate(theta, parmlist[[i]])
         if (!inherits(temp, 'bdsmatrix')) 
             stop("penalty matrix for term", i, "must be a bdsmatrix")
         if (any(dim(temp)) != rep(fcount[i] + zcount[i],2))
             stop ("Invalid penalty matrix for term", i)
         if (fcount2[i] >0){
             t1 <- 1:fcount2[i]; t2 <- 1:fcount[i]
             R[t2 + indx2[i], t1+indx1[i]] <- temp@rmat[t2, t1]
             }
         if (zcount[i] >0) {
             t1 <-  1:zcount[i];  t2 <- t1 + fcount2[i]
             R[indx4[i] + t1, indx3[i]+t1] <- temp@rmat[t2, t2]
             }
         if (i==1) tsave <- temp
         }
     
     bdsmatrix(blocksize=tsave@blocksize, blocks=tsave@blocks, R=R)
     }    
 logfun <- function(theta, varlist, parmlist, kfun,
                    init, fit0, iter, ofile) {
     gkmat <- gchol(kfun(theta, varlist, parmlist))
     ikmat <- solve(gkmat)  #inverse of kmat, which is the penalty
     if (any(diag(ikmat) <=0)) { #Not an spd matrix
         return(0)  # return a "worse than null" fit
         }
     fit <- .C(ofile,
               iter= as.integer(c(iter,iter)),
               beta = as.double(init),
               loglik = double(2),
               as.double(ikmat@blocks),
               as.double(ikmat@rmat),
               hdet = double(1))
     ilik <- fit$loglik[2] -
              .5*(sum(log(diag(gkmat))) + fit$hdet)

     -(1+ ilik - fit0)
     }
 if (length(theta)) {
     logpar <- list(varlist=varlist, parmlist=parmlist,
                    kfun=kfun, init=c(rep(0., npenal), fit0$coef),
                    fit0= fit0$loglik[2],
                    iter=coxme.control$inner.iter,
                    ofile=ofile)
   
     mfit <- do.call('optim', c(list(par= theta[!fixed], fn=logfun, gr=NULL), 
                            optpar, logpar))
     theta <- mfit$par
     }
     gkmat <- gchol(kfun(theta, varlist, parmlist))
     ikmat <- solve(gkmat)  #inverse of kmat, which is the penalty
     fit <- .C(ofile,
               iter= as.integer(c(0, control$iter.max)),
               beta = as.double(c(rep(0., nfrail), fit0$coef)),
               loglik = double(2),
               as.double(ikmat@blocks),
               as.double(ikmat@rmat),
               hdet = double(1))
     ilik <- fit$loglik[2] -
              .5*(sum(log(diag(gkmat))) + fit$hdet)
  if (refine.n > 0) {
      bmat <- matrix(rnorm(length(beta)*refine.n), ncol=refine.n)
      sim.pen <- colSums(bmat^2)/2   #This is b' \Sigma^{-1} b /2 

      rfit <- .C(rfun,
                 as.integer(refine.n),
                 as.double(beta),
                 as.double(gkmat %*% bmat),
                 loglik = double(refine.n),
                 approx = double(refine.n))
      errhat <- exp(rfit$loglik- ilik) - 
                        exp(fit$loglik[2] + sim.pen - (rfit$approx + ilik))
      ilik = ilik + log(1 + mean(errhat))
      r.correct <- c(correction=mean(errhat), var=var(errhat)/refine.n)
      }
  nfrail <- nrow(ikmat)  #total number of penalized terms
  nvar2  <- nvar + (nfrail - nsparse)  # total number of non-sparse coefs
  nvar3  <- nvar + nfrail              # total number of coefficients
  btot   <- length(ikmat@blocks)

  fit3 <- .C('coxfit6c',
                 u    = double(nvar3),
                 h.b  = double(btot),
                 h.r  = double(nvar2*nvar3),
                 hi.b = double(btot),
                 hi.r = double(nvar2*nvar3),
                 hrank= integer(1),
                 as.integer(ncol(y)),
                 )
      if (nvar2 ==0) {
          hmat <- new('gchol.bdsmatrix', .Dim=c(nvar3, nvar3),
                      blocksize=ikmat@blocksize, blocks=fit3$h.b,
                      rmat=numeric(0), rank=fit3$hrank)
          hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b)
          }
      else {
          rmat1 <- matrix(fit3$h.r, nrow=nvar3)
          rmat2 <- matrix(fit3$hi.r, nrow=nvar3)
          if (nvar ==1 ) {
              rmat1[nvar3,] <- rmat1[nvar3,]/scale
              rmat2[nvar3,] <- rmat2[nvar3,]/scale
              rmat1[,nvar2] <- rmat1[,nvar2]*scale
              rmat2[,nvar2] <- rmat2[,nvar2]/scale
              rmat1[nvar3,nvar2] <- rmat1[nvar3,nvar2]*scale^2
              u <- fit3$u  # the efficient score vector U
              u[nvar3] <- u[nvar3]*scale
              }
          else if (nvar >1) {
              temp <- (nvar3-nvar):nvar3 
              u <- fit$u
              u[temp] <- u[temp]*scale
              rmat1[temp,] <- (1/scale)*rmat1[temp,] #multiply rows* scale
              rmat2[temp,] <- (1/scale)*rmat2[temp,] 

              temp <- (nvar2-nvar):nvar2        #multiply cols
              rmat1[,temp] <- rmat1[,temp] %*% diag(scale)
              rmat2[,temp] <- rmat2[,temp] %*% diag(1/scale)
              temp <- seq(length=length(scale), to=length(rmat1), by=1+nvar3)
              rmat1[temp] <- rmat1[temp]*(scale^2)    #fix the diagonal
              }
          hmat <- new('gchol.bdsmatrix', .Dim=c(nvar3, nvar3),
                      blocksize=ikmat@blocksize, blocks=fit3$h.b,
                      rmat= rmat1, rank=fit3$hrank)
          hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b,
                            rmat=rmat2)
          }
  traceprod <- function(H, P) {
      #block-diagonal portions match
      nfrail <- nrow(P)  #penalty matrix
      temp1 <- sum(H@blocks * P@blocks)
      if (length(P@rmat) >0) {
          #I only want the penalized part of H
          rd <- dim(P@rmat)
          temp1 <- temp1 + sum(H@rmat[1:rd[1], 1:rd[2]] * P@rmat)
          }
      2*temp1 - sum(diag(H)[1:nfrail] * diag(P))
      }
  df <- nvar + (nfrail - traceprod(hinv, ikmat))
      penalty <- sum(fcoef * (ikmat %*% fcoef))
      idf <- nvar + sum(ntheta)

      if (nvar > 0) {
          out <- list(coefficients=list(fixed=fit$beta[-(1:nfrail)]/scale, 
                                 random=theta),
               frail=fit$beta[1:nfrail], penalty=penalty,
               loglik=c(fit0$log[1], ilik, fit$log[2]), var=hinv,
               df=c(idf, df), hmat=hmat, iter=iter, control=control,
               u=u, means=means, scale=scale)
          }
      else out <- list(coefficients=list(fixed=NULL, random=theta),
                frail=fit$beta[1:nfrail], penalty= penalty,
                loglik=c(fit0$log[1], ilik, fit$log[2]), var=hinv,
                df=c(idf, df), hmat=hmat, iter=iter, control=control,
                u=fit3$u, means=means, scale=scale)    

      if (refine.n>0) out<- c(out, list(errhat=errhat, refine=r.correct))

      out
      }
