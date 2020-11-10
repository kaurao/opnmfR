# modified from: https://github.com/asotiras/brainparts

# test run on "iris" data
#' @export
opnmfR_test <- function(r=2, W0="nndsvd", ...) {
  data("iris")
  cat("opnmfR ")
  start_time <- Sys.time()
  nn <- opnmfR(data.matrix(iris[,1:4]), r=r, W0=W0, ...)
  end_time <- Sys.time()
  show(end_time - start_time)
  show(nn$H)
  plot(nn$W[,1], nn$W[,2], col=iris$Species, xlab="Factor 1", ylab="Factor 2")
  
  # now using rcpp
  cat("opnmfRcpp ")
  start_time <- Sys.time()
  nn <- opnmfRcpp(data.matrix(iris[,1:4]), r=r, W0=W0, ...)
  end_time <- Sys.time()
  show(end_time - start_time)
  show(nn$H)
  
  par(mfrow=c(1,2))
  image(nn$W, main="W")
  image(nn$H, main="H")
}

#' @export
opnmfRcpp <- function(X,r,W0=NULL,max.iter=50000,tol=1e-5,memsave=TRUE,eps=1e-16,...) {
  start_time <- Sys.time()
  if(is.null(W0) || is.character(W0)) W0 <- opnmfR_init(X, r, W0)
  W0[W0 < eps] <- 0
  
  res <- opnmfR:::opnmfR_opnmfRcpp(X, W0, tol, max.iter, eps, memsave)
  
  post <- opnmfR_postprocess(X, res$W)
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(W=post$W, H=post$H, diffW=res$diffW, iter=res$iter, mse=post$mse, W0=W0, time=tot_time))
}

#' @export
opnmfR <- function(X,r,W0=NULL,max.iter=50000,tol=1e-5,memsave=TRUE,eps=1e-16, use.gpu=FALSE, ...) {
  # mem=FALSE is faster for large data
  start_time <- Sys.time()
  
  if(any(is.na(X))) {
    Xmean <- median(X, na.rm=TRUE)
    X[is.na(X)] <- Xmean
  }
  
  if(is.null(W0) || is.character(W0)) W0 <- opnmfR_init(X, r, W0)
  W0[W0 < eps] <- 0
  
  W <- W0
  if(!memsave) { XX <- X %*% t(X) }
  

  if(use.gpu) {
     library(gpuR)
     X <- gpuMatrix(X)
     W <- gpuMatrix(W)
     if(!memsave) XX <- gpuMatrix(XX)
  }
 
  pbar <- txtProgressBar(min=0, max=max.iter, style=1)
  for(iter in 1:max.iter) {
    Wold <- W

    # two "types" of update rules
    # "simplified" update rule using XXW
    if(!memsave) XXW <- XX %*% W
    else XXW <- X %*% crossprod(X,W)
    
    W <- W * XXW / (W  %*% crossprod(W,XXW))
    
    # these are the "original" update rules
    #if(!memsave) W <- W * (XX %*% W) / (W %*% (t(W) %*% XX %*% W))
    #else W <- W * (X %*% (t(X) %*% W)) / (W %*% ((t(W) %*% X) %*% (t(X) %*% W))) 
   
    if(!use.gpu) {
      W[W < eps] <- eps
      W <- W / norm(W,"2") # spectral norm
    } else {
      n2 <- svd(t(W) %*% W, 1, 1)
      n2 <- sqrt(as.vector(n2$d)[1])
      W <- W / n2
    }
    
    diffW <- norm(Wold-W, "F") / norm(Wold, "F")
    setTxtProgressBar(pbar, iter)
    if(diffW < tol) break
  }
  close(pbar)
 
  if(use.gpu) W <- as.matrix(W)
 
  post <- opnmfR_postprocess(X, W)
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(W=post$W, H=post$H, diffW=diffW, iter=iter, mse=post$mse, W0=W0, time=tot_time))
}

#' @export
opnmfR_postprocess <- function(X, W) {
  # get H
  H <- t(W) %*% X
  
  # reorder
  hlen <- sqrt(apply(H, 1, function(xx) sum(xx^2)))
  if(any(hlen==0)) cat('WARN low rank:', sum(hlen==0), 'factors have norm 0')
  hlen[hlen==0]  <- 1
  
  Wh <- W
  for(i in 1:ncol(W)) Wh[,i] <- Wh[,i] * hlen[i]
  Whlen <- apply(Wh, 2, function(xx) sum(xx^2))
  ix <- sort(Whlen, decreasing = TRUE, index.return=TRUE)$ix
  W <- W[,ix,drop=FALSE]
  
  # get H from reordered W
  H <- t(W) %*% X
  
  # reconstruction error
  mse <- norm(X-(W %*% H), "F")
  
  return(list(W=W, H=H, mse=mse))
}

#' @export
opnmfR_mse <- function(X, W, H) {
  return(norm(X-(W %*% H), "F"))
}


chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

#' @export
opnmfR_ranksel_ooserr <- function(X, rs, W0=NULL, use.rcpp=TRUE, nrepeat=1, nfold=2, plots=TRUE, seed=NA, ...) {
  stopifnot(ncol(X)>=max(rs))
  start_time <- Sys.time()
  
  if(is.na(seed)) seed <- sample(1:10^6, 1)
  mse <- list()
  nrun <- 0
  for(n in 1:nrepeat) {
    mse[[n]] <- list()
    mse[[n]]$train <- matrix(NA, nfold, length(rs))
    mse[[n]]$test <- matrix(NA, nfold, length(rs))
    mse[[n]]$test_delta <- matrix(NA, nfold, length(rs))
    mse[[n]]$test_test <- matrix(NA, nfold, length(rs))
    mse[[n]]$test_train <- matrix(NA, nfold, length(rs))
    # get folds
    set.seed(seed+n)
    folds <- chunk(sample(1:nrow(X)), nfold)
    for(f in 1:nfold) {
      testidx <- folds[[f]]
      trainidx <- setdiff(1:nrow(X), testidx)
      cat("repeat", n, ": fold", f, ": nex",  length(testidx), ":")
      
      factrain <- list()
      factest <- list()
      for(r in 1:length(rs)) {
        if(use.rcpp) {
          factrain[[r]] <- opnmfRcpp(X[trainidx,], rs[r], W0=W0, ...)
          factest[[r]] <- opnmfRcpp(X[testidx,], rs[r], W0=W0, ...)
        } else {
          factrain[[r]] <- opnmfR(X[trainidx,], rs[r], W0=W0, ...)
          factest[[r]] <- opnmfR(X[testidx,], rs[r], W0=W0, ...)
        }
        
        mse[[n]]$train[f,r] <- factrain[[r]]$mse
        
        # project on train H and get reconstruction error
        #Xtest_test <- factest[[r]]$W %*% factest[[r]]$H
        #Xtest_train <- factest[[r]]$W %*% factrain[[r]]$H
        xtest_test <- opnmfR_reconstruct(X[testidx,], factest[[r]]$H)
        xtest_train <- opnmfR_reconstruct(X[testidx,], factrain[[r]]$H)
        mse[[n]]$test[f,r] <- norm(xtest_test - xtest_train, "F")
        mse[[n]]$test_test[f,r] <- norm(xtest_test - X[testidx,], "F")
        mse[[n]]$test_train[f,r] <- norm(xtest_train - X[testidx,], "F")
        mse[[n]]$test_delta[f,r] <- abs(norm(X[testidx,] - xtest_train, "F") - norm(xtest_test - X[testidx,], "F"))
        
        nrun <- nrun + 1
        cat("test err", mse[[n]]$test[f,r], "\n#######\n")
      } # rs
    } # nfold
  } # nrepeat
  
  errtrain <- do.call(rbind, lapply(mse, function(z) apply(z$train, 2, mean)))
  errtest <- do.call(rbind, lapply(mse, function(z) apply(z$test, 2, mean)))
  errtest_delta <- do.call(rbind, lapply(mse, function(z) apply(z$test_delta, 2, mean)))
  
  mse$train <- errtrain
  mse$test <- errtest
  mse$test_delta <- errtest_delta
  rownames(mse$train) <- paste("repeat", 1:nrepeat, sep="")
  colnames(mse$train) <- paste("rank", rs, sep="")
  rownames(mse$test) <- paste("repeat", 1:nrepeat, sep="")
  colnames(mse$test) <- paste("rank", rs, sep="")
  rownames(mse$test_delta) <- paste("repeat", 1:nrepeat, sep="")
  colnames(mse$test_delta) <- paste("rank", rs, sep="")
  
  errtrain <- apply(errtrain, 2, mean)
  errtest <- apply(errtest, 2, mean)
  errtest_delta <- apply(errtest_delta, 2, mean)
  selr <- rs[which.min(errtest)]
  selr_delta <- rs[which.min(errtest_delta)]
  
  if(plots) {
    par(mfrow=c(1,3))
    plot(errtrain, main="Train", ylab="MSE", xlab="Rank")
    plot(errtest, main="Test", ylab="MSE", xlab="Rank")
    points(selr, errtest[selr], cex=1.5, col="red", pch=10)
    plot(errtest_delta, main="Test", ylab="MSE (delta)", xlab="Rank")
    points(selr_delta, errtest_delta[selr_delta], cex=1.5, col="red", pch=10)
  }
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(mse=mse, time=tot_time, seed=seed, selected=selr, selected_delta=selr_delta))
}


#' @export
opnmfR_reconstruct <- function(X, H) {
  library(MASS) # for ginv
  X <- abs(X %*% (t(H) %*% ginv(H %*% t(H))) %*% H)
  return(X)
}


#' @export
opnmfR_ranksel_perm <- function(X, rs, W0=NULL, use.rcpp=TRUE, nperm=1, plots=TRUE, seed=NA, ...) {
  stopifnot(ncol(X)>=max(rs))
  start_time <- Sys.time()
  
  # original data
  cat("original data... ")
  nn <- list()
  mse <- list()
  for(r in 1:length(rs)) {
    if(use.rcpp) {
      nn[[r]] <- opnmfRcpp(X, rs[r], W0=W0, ...)
    } else {
      nn[[r]] <- opnmfR(X, rs[r], W0=W0, ...)
    }
    
    cat("orig", "rank:", r, "iter:", nn[[r]]$iter, "time:", nn[[r]]$time, "diffW:", nn[[r]]$diffW, "\n")

    mse[[r]] <- list()
    mse[[r]]$orig <- opnmfR_mse(X, nn[[r]]$W, nn[[r]]$H)
  }
  cat("done\n")
  
  # permuted data
  cat("permuted data... ")
  if(is.na(seed)) seed <- sample(1:10^6, 1)
  for(p in 1:nperm) {
    set.seed(seed+p)
    Xperm <- apply(X,2,sample) # permute the rows in each column
    for(r in 1:length(rs)) {
      if(use.rcpp) {
        nnp <- opnmfRcpp(Xperm, rs[r], W0=W0, ...)
      } else {
        nnp <- opnmfR(Xperm, rs[r], W0=W0, ...)
      }
    
      cat("perm", p, "rank:", r, "iter:", nnp$iter, "time:", nnp$time,"diffW:", nnp$diffW, "\n")
      mse[[r]]$perm <- c(mse[[r]]$perm, opnmfR_mse(Xperm, nnp$W, nnp$H))
    }
  }
  cat("done\n")
  
  names(mse) <- rs
  names(nn) <- rs
  
  mseorig <- sapply(mse, function(xx) xx$orig)
  mseperm <- sapply(mse, function(xx) mean(xx$perm))
  mseorig <- mseorig / max(mseorig)
  mseperm <- mseperm / max(mseperm)
  
  sel <- which(diff(mseorig) > diff(mseperm))
  selr <- rs[sel]
  
  if(plots) {
    maxi <- max(c(mseorig, mseperm))
    mini <- min(c(mseorig, mseperm))
    
    plot(NA, xlim=c(min(rs),max(rs)), ylim=c(mini,maxi), xlab="Rank", ylab="MSE (scaled)")
    points(rs, mseorig, type='b', pch=16)
    points(rs, mseperm, type='b', pch=17, lty=2)
    points(selr, mseorig[sel], cex=2, pch=1, lwd=2, col="red")
    legend("topright", legend = c("Orig.","Perm."), pch=16:17)
    title(main="Permutation based rank selection", sub=paste("Selected rank = ", min(selr)))
  }
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(mse=mse, selected=selr, factorization=nn[sel], time=tot_time, seed=seed))
  
}

#' @export
opnmfR_init <- function(X, r, W0) {
  cat("initializing... r =", r, " ")
  n <- nrow(X)
  if(is.null(W0) || W0=="random") {
    cat("random")
    W0 <- matrix(runif(n*r), ncol = r)
  } else if(is.character(W0) && W0=="nndsvd") {
    cat(W0)
    W0 <- opnmfR_nndsvd(X,r)$W
  } else if(is.character(W0) && W0=="nndsvd1") {
    cat(W0)
    W0 <- opnmfR_nndsvd(X,r,flag = 1)$W
  } else if(is.character(W0) && W0=="nndsvd2") {
    cat(W0)
    W0 <- opnmfR_nndsvd(X,r,flag = 2)$W
  } else if(is.character(W0) && W0=="brainparts") {
    cat(W0)
    W0 <- opnmfR_nndsvd(X,r,epseps=0.1)$W
  }
  cat(" done\n")
  return(W0)
}


# based on matlab code from: http://www.boutsidis.org/software.html
# to make it similar to brainparts set epseps=0.1
#' @export
opnmfR_nndsvd  <- function(x,k,flag=0,seed=NA,s=NULL,eps=0.0000000001,epseps=0) {
  x[is.na(x)] <- 0
  stopifnot(all(x>=0))
  
  if(is.null(s)) s <- svd(x,k,k)
  stopifnot(ncol(s$u)>=k)
  stopifnot(ncol(s$v)>=k)
  
  m <- nrow(x)
  n <- ncol(x)
  W <- matrix(0,nrow=m,ncol=k)
  H <- matrix(0,nrow=k,ncol=n)
  W[,1] <- sqrt(s$d[1])*abs(s$u[,1])
  H[1,] <- sqrt(s$d[1])*abs(s$v[,1])
  
  if(k>1) {
    for(i in seq(2,k)) {
      uu <- s$u[,i]; vv <- s$v[,i]
      uup <- (uu>=0) * uu; uun <- (uu<0) * (-uu)
      vvp <- (vv>=0) * vv; vvn <- (vv<0) * (-vv)
      n_uup <- sqrt(sum(uup^2)); n_uun <- sqrt(sum(uun^2))
      n_vvp <- sqrt(sum(vvp^2)); n_vvn <- sqrt(sum(vvn^2))
      termp <- n_uup*n_vvp; termn <- n_uun*n_vvn
      if(termp >= termn) {
        W[,i] <- sqrt(s$d[i]*termp)*uup/n_uup
        H[i,] <- sqrt(s$d[i]*termp)*vvp/n_vvp
      } else {
        W[,i] <- sqrt(s$d[i]*termn)*uun/n_uun
        H[i,] <- sqrt(s$d[i]*termn)*vvn/n_vvn
      }
    }
  }
  
  # these should be set to 0 but brainparts sets it to 0.1
  W[W<=eps] <- epseps
  H[H<=eps] <- epseps
  
  # these flags have no effect as above we set small numbers to 0.1
  if(flag==1) {
    average <- mean(x)
    W[W==0] <- average
    H[H==0] <- average
  } else if(flag==2) {
    if(is.na(seed)) seed <- sample(1:10^6,1)
    set.seed(seed)
    average <- mean(x)
    ind1 <- W==0
    ind2 <- H==0
    W[ind1] <- average*runif(sum(ind1))/100
    H[ind2] <- average*runif(sum(ind2))/100
  } else if(flag>0) stop("unknown flag in opnmfR_nndsvd")
  
  return(list(W=W, H=H, seed=seed))
  
}

