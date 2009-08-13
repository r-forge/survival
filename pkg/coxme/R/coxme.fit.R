# Automatically generated from all.nw using noweb
 coxme.fit <- function(x, y, strata, offset, ifixed, control,
                         weights, ties, rownames, 
                         fmat, zmat, varlist, vparm, theta,
                         ntheta, refine.n) {
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
 if (is.null(ifixed) ) ifixed <- rep(0., ncol(x))
 else if (length(ifixed) != ncol(x))
     stop("Wrong lenth for initial parameters of the fixed effects")
                                   
 fit0 <- coxph(y ~ x + offset(offset), weights=weights, method=ties,
               init=ifixed, iter=0)
 kfun <- function(theta, varlist, vparm, ntheta,
                  fcount=apply(fmat,2,max), zcount = ncol(Z)) {

     nrandom <- length(varlist)
     sindex <- rep(1:nrandom, ntheta) #which thetas to which terms

     tmat <- varlist[[1]]$generate(theta[sindex==1], vparm[[1]]) 
     if (!inherits(tmat, 'bdsmatrix')) 
         tmat <- bdsmatrix(blocksize=integer(0), blocks=numeric(0), rmat=tmat)
         
     if (nrandom ==1) return(tmat)

     # Need to build up the matrix by pasting up a composite R
     nsparse <- sum(tmat@blocksize)
     nrow.R <- sum(fcount) + zcount
     ncol.R <- nrow.R - nsparse
     R <- matrix(0., nrow.R, ncol.R)
     indx1 <- 0            #current row  wrt filling in intercepts
     indx2 <- sum(fcount)  #current row wrt filling in slopes
     
     if (ncol(tmat) > nsparse) {
         k <- (nsparse+1):ncol(tmat)
         temp <- as.matrix(tmat[k,k])
         }
     if (fcount[1] > nsparse) {
         j <- fcount[1] - nsparse   #number of intercept columns
         R[1:nrow(temp), 1:j] <- temp[,1:j]
         indx1 <- indx1 +j
         }
     else j <- 0
     if (ncol(tmat) > fcount[1]) {
         k <- ncol(tmat) - fcount[1]  #number of slopes
         if (fcount[1]>0) {
             R[1:fcount[1], indx2+1:k] <- temp[1:fcount[1], 1:k]
             R[1:k, indx2+1:k] <- temp[fcount[1]+1:k, 1:k]
             }
         else R[1:nrow(temp1) + indx2, 1:k+indx2] <- temp[, j +1:k]

         indx2 <- indx2+k
         }
     
     for (i in 2:nrandom) {
         temp <- as.matrix(varlist[[i]]$generate(theta[sindex==i], vparm[[i]]))
         j <- fcount[i]
         if (j >0) {
             R[1:j + indx1, 1:j+indx1] <- temp[1:j, 1:j]
             }
         if (ncol(temp) > fcount[i]) {
             k <- ncol(temp) - fcount[i]
             if (j>0) 
                 R[1:j+indx1, 1:k + indx2] <- temp[1:j, j+1:k]
             R[1:k+indx2, 1:k + indx2] <- temp[j+ 1:k, j+ 1:k]
             indx2 <- indx2 +k
             }
         indx1 <- indx1 +j
         }
     
     bdsmatrix(blocksize=tmat@blocksize, blocks=tmat@blocks, R=R)
     }    
 dummy <- kfun(theta, varlist, vparm, ntheta)
 if (is.null(dummy@rmat)) rcol <- 0
     else                 rcol <- ncol(dummy@rmat)
 npenal <- ncol(dummy)  #total number of penalized terms

 if (ncol(fmat) >1) {
     for (i in 2:ncol(fmat)) fmat[,i] <- fmat[,i] + max(fmat[,i-1])
     }  
 if (ncol(fmat)>0) {
     findex <- matrix(0, nrow=max(fmat), ncol=ncol(fmat))
     for (i in 1:ncol(fmat)) findex[cbind(fmat[,i], i)] <- 1
     }
 else findex <- 0  # dummy value
     
 if (is.null(control$sparse.calc)) {
     nevent <- sum(y[,ncol(y)])
     if (length(dummy@blocksize)<=1) nsparse<- 0
     else nsparse <- sum(dummy@blocksize)
     itemp <- max(c(0,fmat)) - nsparse  #number of non-sparse intercepts
     
     if ((2*n) > (nevent*(nsparse-itemp))) control$sparse.calc <- 0
     else control$sparse.calc <- 1
     }

 ifit <- .C('coxfit6a', 
                as.integer(n),
                as.integer(nvar),
                as.integer(ncol(y)),
                as.double(c(y)),
                as.double(cbind(zmat,x)),
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
                as.integer(rcol),
                means = double(nvar),
                scale = double(nvar),
                as.integer(ties=='efron'),
                as.double(control$toler.chol),
                as.double(control$eps),
                as.integer(control$sparse.calc))
     means   <- ifit$means
     scale   <- ifit$scale
 logfun <- function(theta, varlist, vparm, kfun, ntheta,
                    init, fit0, iter, ofile) {
     gkmat <- gchol(kfun(theta, varlist, vparm, ntheta))
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
     logpar <- list(varlist=varlist, vparm=vparm, 
                    ntheta=ntheta, kfun=kfun,
                    init=c(rep(0., npenal), scale*fit0$coef),
                    fit0= fit0$loglik[2],
                    iter=control$inner.iter,
                    ofile=ofile)
   
     mfit <- do.call('optim', c(list(par= theta, fn=logfun, gr=NULL), 
                            control$optpar, logpar))
     theta <- mfit$par
     iter <- mfit$counts[1] * c(1, control$inner.iter)
     }
 else iter <- c(0,0)

     gkmat <- gchol(kfun(theta, varlist, vparm, ntheta))
     ikmat <- solve(gkmat)  #inverse of kmat, which is the penalty
     fit <- .C(ofile,
               iter= as.integer(c(0, control$iter.max)),
               beta = as.double(c(rep(0., npenal), fit0$coef*scale)),
               loglik = double(2),
               as.double(ikmat@blocks),
               as.double(c(ikmat@rmat,0)),
               hdet = double(1))
     ilik <- fit$loglik[2] -
              .5*(sum(log(diag(gkmat))) + fit$hdet)
     iter[2] <- iter[2] + fit$iter[2]
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
  nsparse <- sum(ikmat@blocksize)
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
                 as.integer(ncol(y))
                 )
      if (nvar2 ==0) {
          hmat <- new('gchol.bdsmatrix', Dim=c(nvar3, nvar3),
                      blocksize=ikmat@blocksize, blocks=fit3$h.b,
                      rmat=matrix(0,0,0), rank=fit3$hrank,
                      Dimnames=list(NULL, NULL))
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
              temp <- seq(to=nvar3, length=length(scale))
              u <- fit3$u
              u[temp] <- u[temp]*scale
              rmat1[temp,] <- (1/scale)*rmat1[temp,] #multiply rows* scale
              rmat2[temp,] <- (1/scale)*rmat2[temp,] 

              temp <- temp-nsparse          #multiply cols
              rmat1[,temp] <- rmat1[,temp] %*% diag(scale)
              rmat2[,temp] <- rmat2[,temp] %*% diag(1/scale)
              temp <- seq(length=length(scale), to=length(rmat1), by=1+nvar3)
              rmat1[temp] <- rmat1[temp]*(scale^2)    #fix the diagonal
              }
          hmat <- new('gchol.bdsmatrix', Dim=c(nvar3, nvar3),
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
  newtheta <- NULL   
  nrandom <- length(varlist)
  sindex <- rep(1:nrandom, ntheta) #which thetas to which terms

  for (i in 1:nrandom) 
      newtheta <- c(newtheta, varlist[[i]]$wrapup(theta[sindex==i], vparm[[i]]))
      fcoef <- fit$beta[1:nfrail]
      penalty <- sum(fcoef * (ikmat %*% fcoef))/2
      idf <- nvar + sum(ntheta)

      if (nvar > 0) {
          out <- list(coefficients=list(fixed=fit$beta[-(1:nfrail)]/scale, 
                                 random=newtheta),
               frail=fit$beta[1:nfrail], penalty=penalty,
               loglik=c(fit0$log[1], ilik, fit$log[2]), var=hinv,
               df=c(idf, df), hmat=hmat, iter=iter, control=control,
               u=u, means=means, scale=scale)
          }
      else out <- list(coefficients=list(fixed=NULL, random=newtheta),
                frail=fit$beta[1:nfrail], penalty= penalty,
                loglik=c(fit0$log[1], ilik, fit$log[2]), var=hinv,
                df=c(idf, df), hmat=hmat, iter=iter, control=control,
                u=fit3$u, means=means, scale=scale)    

      if (refine.n>0) out<- c(out, list(errhat=errhat, refine=r.correct))

      out
      }
