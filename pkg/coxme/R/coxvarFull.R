# Automatically generated from all.nw using noweb

coxvarFull <- function(collapse=TRUE) {
    collapse <- collapse
    # Because of R's lexical scoping, the values of the options
    #  above, at the time the line below is run, are known to the
    #  init function
    init <- function(init, fixed, Z, G, gnames, sparse) {
        ncluster <- length(gnames) # also = ncol(G)
        my.init <- rep(.2, ncluster)
        names(my.init) <- gnames
        my.fixed <- rep(0, ncluster)

        if (!missing(init)) {
            if (length(init) != ncluster) 
                    stop ("Wrong length for initial values")
            my.init<- init
            }

        if (!missing(fixed)) {
            if (length(fixed) != ncluster)
                stop ("Wrong length for fixed variance vector")
            my.fixed <- fixed
            }
        if (is.null(Z)) nx <-0 else nx <- ncol(Z))
        if (ncluster==1) {
            temp <- table(G)
            ngroup <- length(temp) * (nx+1)
            nx2 <- nx* ngroup
            vmat <
            if (length(temp) >= sparse[1] && any(temp/sum(temp) < sparse[2])){
                which.sparse <- which(temp/sum(temp) < sparse[2])
                n.sparse <- length(which.sparse)
                index <- c(which.sparse, 1:ngroup[-which.sparse])
                n.dense <- ngroup - n.sparse
                if (n.dense > 1) {
                    vmat <- bdsmatrix(blocksize=rep(1, n.sparse),
                                      blocks=rep(1.0, nsparse),
                                  rmat=rbind(matrix(0., n.sparse, n.dense+nx2),
                                             diag(n.dense+ nx2)))
                    }
                else {
                    if (nx==0) vmat <- bdsmatrix(blocksize=rep(1, ngroup),
                                                 blocks = rep(1.0, ngroup))
                    else vmat <- bdsmatrix(blocksize=rep(1, ngroup),
                                           blocks = rep(1.0, ngroup),
                                           rmat= rbind(matrix(0., ngroup,nx2),
                                                       diag(nx2))
                }
            else {
                index <- 1:ngroup
                vmat <- bdsmatrix(blocksize=ngroup, blocks=diag(ngroup))
                }
          list(theta=my.init,  
                

                
        varmat <- vector(ncluster, 'list') 
        if (ncluster >1) {
            if (collapse) {
                g <- strata(G, shortlabel=TRUE)
                
                
            temp <- groups
            for (i in 2:ncluster)
                    temp[,i] <- strata2(groups[,1:i], shortlabel=shortlabel,
                                 sep='/')
            groups <- temp
            }
