# partly based on `brainparts`: https://github.com/asotiras/brainparts

#' opnmfR: A package for orthonormal projective non-negative matrix factorization
#'
#' This package provides functions for opnmf factorization using R code as well as Rcpp code as well as 
#' for rank selection and initializatoin.
#' @docType package
#' @name opnmfR
#' @useDynLib opnmfR, .registration=TRUE
NULL
#> NULL


library(NMF)
library(lpSolve)
library(aricode)

#' Simple test
#'
#' @param X A matrix, if NULL the "iris" data is used (default NULL)
#' @param r A number, rank to use (default 2)
#' @param W0 A string or matrix for initialization (default "nndsvd")
#' @return A list with factorization results using R and Rcpp function calls
#' @examples
#' result <- opnmfR_test()
#' @export
opnmfR_test <- function(X=NULL, y=NULL, r=2, W0="nndsvd", ...) {
  if(is.null(X)) {
    data("iris")
    X <- data.matrix(iris[,1:4])
    y <- iris$Species
  }
  
  # using R
  cat("opnmfR test ")
  start_time <- Sys.time()
  nn <- opnmfR(X, r=r, W0=W0, ...)
  end_time <- Sys.time()
  show(end_time - start_time)
  show(nn$H)
  if(r>=2) {
    plot(nn$W[,1], nn$W[,2], col=y, xlab="Factor 1", ylab="Factor 2")
  }
  
  # using rcpp
  cat("opnmfRcpp test ")
  start_time <- Sys.time()
  nncpp <- opnmfRcpp(X, r=r, W0=W0, ...)
  end_time <- Sys.time()
  show(end_time - start_time)
  show(nn$H)
  
  par(mfrow=c(1,2))
  image(nn$W, main="W")
  image(nn$H, main="H")
  par(mfrow=c(1,1))
  
  return(list(nn=nn, nncpp=nncpp))
}

#' @export
opnmfR_test_ranksel <- function(X=NULL, rs=NULL, W0=NULL, nrepeat=1) {
  if(is.null(X)) {
    data("iris")
    cat("opnmfR ranksel ")
    start_time <- Sys.time()
    X <- data.matrix(iris[,1:4])
    X <- t(cbind(X, X)) # duplicate columns and transpose, i.e. factorize features
  }
  
  if(is.null(rs)) rs <- 1:nrow(X)
  
  perm <- opnmfR_ranksel_perm(X, rs, W0=W0, nperm=nrepeat)
  ooser <- opnmfR_ranksel_ooser(X, rs, W0=W0)
  splithalf <- opnmfR_ranksel_splithalf(X, rs, W0=W0, nrepeat=nrepeat)
  
  return(list(perm=perm, ooser=ooser, splithalf=splithalf))
}

#' @export
opnmfR_test_ranksel_synthetic <- function(n=100, r=10, p=100, rs=NULL, W0=NULL, nrepeat=1) {
  stopifnot(r <= 20) # keep the test computable
  stopifnot(r < n && r < p)
  X <- syntheticNMF(n,r,p)
  if(is.null(rs)) rs <- 1:min(nrow(X), round(r+r*0.5))

  perm <- opnmfR_ranksel_perm(X, rs, W0=W0, nperm=nrepeat, rtrue=r)
  ooser <- opnmfR_ranksel_ooser(X, rs, W0=W0, rtrue=r)
  splithalf <- opnmfR_ranksel_splithalf(X, rs, W0=W0, nrepeat=nrepeat, rtrue=r)
  
  return(list(perm=perm, ooser=ooser, splithalf=splithalf))
}

#' @export
opnmfR_compare_nmf <- function(n=100, r=10, p=100, nmfalgo="snmf/r", rs=NA) {
  stopifnot(r < n && r < p)
  if(is.na(rs)) rs <- r
  stopifnot(rs < n && rs < p)
  X <- syntheticNMF(n,r,p)
  sim_cosine <- rep(NA, r)
  ari <- rep(NA, r)
  for(rr in 1:rs) {
    fac1 = nmf(X, rr, nmfalgo)
    fac2 = opnmfR(X, rr)
    sim <- opnmfR_cosine_similarity(t(fac1@fit@W), t(fac2$W))
    lp <- lp.assign(sim, direction = "max")
    sim_cosine[rr] <- lp$objval/rr
    if(rr > 1) {
      ari[rr] <- ARI(apply(fac1@fit@W,1,which.max), apply(fac2$W,1,which.max))
    }
  }
  par(mfrow=c(1,2))
  plot(sim_cosine, type="b", xlab="Rank", ylab="Similarity (cosine)", main=paste("NMF", nmfalgo), pch=16)
  abline(v=r, lty=2, col="gray")
  plot(ari, type="b", xlab="Rank", ylab="Similarity (ARI)", main=paste("NMF", nmfalgo), pch=17)
  abline(v=r, lty=2, col="gray")
  par(mfrow=c(1,1))
  
  return(list(sim_cosine=sim_cosine, ari=ari))
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
opnmfR_ranksel_ooser <- function(X, rs, W0=NULL, use.rcpp=TRUE, nrepeat=1, nfold=2, plots=TRUE, seed=NA, rtrue=NA, ...) {
  # we create folds over the columns
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
    # we need to fold over the columns
    set.seed(seed+n)
    folds <- chunk(sample(1:ncol(X)), nfold)
    for(f in 1:nfold) {
      testidx <- folds[[f]]
      trainidx <- setdiff(1:ncol(X), testidx)
      cat("repeat", n, ": fold", f, ": #test",  length(testidx), ": #train",  length(trainidx), ":")
      
      factrain <- list()
      factest <- list()
      for(r in 1:length(rs)) {
        # get factors
        if(use.rcpp) {
          factrain[[r]] <- opnmfRcpp(X[,trainidx], rs[r], W0=W0, ...)
          factest[[r]] <- opnmfRcpp(X[,testidx], rs[r], W0=W0, ...)
        } else {
          factrain[[r]] <- opnmfR(X[,trainidx], rs[r], W0=W0, ...)
          factest[[r]] <- opnmfR(X[,testidx], rs[r], W0=W0, ...)
        }
        
        # training error
        mse[[n]]$train[f,r] <- factrain[[r]]$mse
        
        # project and get reconstruction error
        #Xtest_test <- factest[[r]]$W %*% factest[[r]]$H
        #Xtest_train <- factest[[r]]$W %*% factrain[[r]]$H
        # reconstruct test data using test data factors
        xtest_test <- opnmfR_reconstruct_W(X[,testidx], factest[[r]]$W)
        # reconstruct test data using training data factors
        xtest_train <- opnmfR_reconstruct_W(X[,testidx], factrain[[r]]$W)
        mse[[n]]$test[f,r] <- norm(xtest_test - xtest_train, "F")
        mse[[n]]$test_test[f,r] <- norm(xtest_test - X[,testidx], "F")
        mse[[n]]$test_train[f,r] <- norm(xtest_train - X[,testidx], "F")
        mse[[n]]$test_delta[f,r] <- abs(norm(X[,testidx] - xtest_train, "F") - norm(xtest_test - X[,testidx], "F"))
        
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
  
  # select rank(s)
  selr_ind <- which.min(errtest)
  selr <- rs[selr_ind]
  selr_delta_ind <- rs[which.min(errtest_delta)]
  selr_delta <- rs[selr_delta_ind]
  
  if(plots) {
    par(mfrow=c(1,3))
    plot(rs, errtrain, main="Train", ylab="MSE", xlab="Rank")
    if(!is.na(rtrue)) abline(v=rtrue, lty=2, col="gray")
    plot(rs, errtest, main="Test", ylab="MSE", xlab="Rank")
    points(selr, errtest[selr_ind], cex=1.5, col="red", pch=10)
    if(!is.na(rtrue)) abline(v=rtrue, lty=2, col="gray")
    plot(rs, errtest_delta, main="Test", ylab="MSE (delta)", xlab="Rank")
    points(selr_delta, errtest_delta[selr_delta_ind], cex=1.5, col="red", pch=10)
    if(!is.na(rtrue)) abline(v=rtrue, lty=2, col="gray")
    par(mfrow=c(1,1))
  }
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(mse=mse, time=tot_time, seed=seed, selected=selr, selected_delta=selr_delta))
}


#' @export
opnmfR_reconstruct_H <- function(X, H) {
  library(MASS) # for ginv
  X <- abs(X %*% (t(H) %*% ginv(H %*% t(H))) %*% H)
  return(X)
}


#' @export
opnmfR_reconstruct_W <- function(X, W) {
  library(MASS) # for ginv
  X <- abs((W %*% ginv(W)) %*% X)
  return(X)
}

opnmfR_cosine_similarity <- function(x, y){
  xy <- tcrossprod(x, y)
  x <- sqrt(apply(x, 1, crossprod))
  y <- sqrt(apply(y, 1, crossprod))
  xy / outer(x,y)
}

#' @export
opnmfR_ranksel_splithalf <- function(X, rs, W0=NULL, use.rcpp=TRUE, nrepeat=1, similarity="inner", splits=NA, plots=TRUE, seed=NA, rtrue=NA, ...) {
  # we create folds over the columns
  library(lpSolve) # for solving the assignment problem
  library(aricode) # for ARI
  
  stopifnot(ncol(X)>=max(rs))
  start_time <- Sys.time()
  
  if(!is.na(splits)) {
    stopifnot(typeof(splits) == "list")
    nrepeat <- length(splits)
  }
  stopifnot(nrepeat > 0)
  
  if(is.na(seed)) seed <- sample(1:10^6, 1)
  perf <- list()
  nrun <- 0
  for(n in 1:nrepeat) {
    perf[[n]] <- list()
    perf[[n]]$train_err <- matrix(NA, 2, length(rs))
    perf[[n]]$sim_inner <- rep(NA, length(rs))
    perf[[n]]$sim_cosine <- rep(NA, length(rs))
    perf[[n]]$ari <- rep(NA, length(rs))
    # cor_cosine is from https://www.sciencedirect.com/science/article/pii/S1053811919309395
    perf[[n]]$cor_cosine <- rep(NA, length(rs))
    
    # get folds
    # we need to fold over the columns
    set.seed(seed+n)
    if(!is.na(splits)) {
      idx1 <- splits[[n]]
    } else {
      folds <- chunk(sample(1:ncol(X)), 2)
      idx1 <- folds[[1]]
    }
    idx2 <- setdiff(1:ncol(X), idx1)
    cat("repeat", n, ": #idx1",  length(idx1), ": #idx2",  length(idx2), ":")
    
    fac1 <- list()
    fac2 <- list()
    for(r in 1:length(rs)) {
      # get factors
      if(use.rcpp) {
        fac1[[r]] <- opnmfRcpp(X[,idx1], rs[r], W0=W0, ...)
        fac2[[r]] <- opnmfRcpp(X[,idx2], rs[r], W0=W0, ...)
      } else {
        fac1[[r]] <- opnmfR(X[,idx1], rs[r], W0=W0, ...)
        fac2[[r]] <- opnmfR(X[,idx2], rs[r], W0=W0, ...)
      }
      
      # training error
      perf[[n]]$train_err[1,r] <- fac1[[r]]$mse
      perf[[n]]$train_err[2,r] <- fac2[[r]]$mse
      
      # match two solutions on their W
      # inner product
      sim <- t(fac1[[r]]$W) %*% fac2[[r]]$W
      lp <- lp.assign(sim, direction = "max")
      perf[[n]]$sim_inner[r] <- lp$objval/rs[r]

      # cosine
      sim <- opnmfR_cosine_similarity(t(fac1[[r]]$W), t(fac2[[r]]$W))
      lp <- lp.assign(sim, direction = "max")
      perf[[n]]$sim_cosine[r] <- lp$objval/rs[r]
      
      if(r > 1) {
        # cor_cosine
        sim1 <- opnmfR_cosine_similarity(fac1[[r]]$W, fac1[[r]]$W) 
        sim2 <- opnmfR_cosine_similarity(fac2[[r]]$W, fac2[[r]]$W)
        stopifnot(nrow(sim1) == nrow(fac1[[r]]$W))
        stopifnot(ncol(sim1) == nrow(fac1[[r]]$W))
        cr <- rep(NA, nrow(sim1))
        for(ii in 1:nrow(sim1)) cr[ii] <- cor(sim1[ii,], sim2[ii,])
        perf[[n]]$cor_cosine[r] <- mean(cr, na.rm = TRUE)
        
        # ari
        cl1 <- apply(fac1[[r]]$W, 1, which.max)
        cl2 <- apply(fac2[[r]]$W, 1, which.max)
        stopifnot(length(cl1)==nrow(X))
        perf[[n]]$ari[r] <- ARI(cl1, cl2)
      }
      
      cat("match", perf[[n]]$sim_cosine[r], "\n#######\n")
    } # rs
  } # nrepeat
  
  selection <- opnmfR_ranksel_splithalf_select(perf, rs, plots, rtrue)
  
  end_time <- Sys.time()
  tot_time <- end_time - start_time
  
  return(list(perf=perf, time=tot_time, seed=seed, selection=selection))
}

opnmfR_ranksel_splithalf_select <- function(perf, rs, plots=TRUE, rtrue=NA) {
  
  measures <- names(perf[[1]])
  measures_avg <- matrix(NA, nrow=length(measures), ncol=length(rs))
  rownames(measures_avg) <- measures
  colnames(measures_avg) <- rs
  
  selected <- matrix(NA, nrow = length(measures), ncol=3)
  rownames(selected) <- measures
  colnames(selected) <- c("value", "index", "rank")
  
  for(i in 1:length(measures)) {
    findmax <- TRUE
    if(measures[i]=="train_err") {
      mea <- do.call(rbind, lapply(perf, function(z) apply(z[[measures[i]]], 2, mean)))
      findmax <- FALSE
    } else{
      mea <- do.call(rbind, lapply(perf, function(z) z[[measures[i]]]))
    }
    
    mea <- apply(mea, 2, mean)
    measures_avg[i,] <- mea
  
    # select rank(s)
    if(findmax) {
      selected[i,1] <- max(mea, na.rm=TRUE)
      selected[i,2] <- which.max(mea==selected[i,1])[1]
      selected[i,3] <- rs[selected[i,2]]
    } else {
      selected[i,1] <- min(mea)
      selected[i,2] <- which.min(mea)
      selected[i,3] <- rs[selected[i,2]]
    }
  }
  
  if(plots) {
    par(mfrow=c(1,1))
    startoverlap <- FALSE
    measures_overlap <- c()
    col_overlap <- c()
    for(i in 1:length(measures)) {
      mea <- measures_avg[i,]
      if(measures[i]=="train_err") {
        plot(mea, xlab="Ranks", ylab="Training MSE", main="Split-half", type="b")
      } else {
        mea <- (mea - min(mea, na.rm=TRUE)) / (max(mea, na.rm=TRUE) - min(mea, na.rm=TRUE))
        if(!startoverlap) {
          par(mar=c(5.3, 4.3, 4.3, 8.3), xpd=TRUE)
          plot(x=rs, y=mea, ylim=c(0,1.1), type="b", col=i, xlab="Ranks", ylab="Measure", main="Split-half")
        }
        else points(x=rs, y=mea, ylim=c(0,1.3), type="b", col=i)
        points(selected[i,3], mea[selected[i,2]], cex=1.5+i/10, col=i, pch=10)
        measures_overlap <- c(measures_overlap, measures[i])
        col_overlap <- c(col_overlap, i)
        startoverlap <- TRUE
      }
    }
    if(!is.na(rtrue)) abline(v=rtrue, lty=2, col="gray")
    legend("topright", inset=c(-0.4,0), legend=measures_overlap, col=col_overlap, pch=1)
  }
  
  return(list(measures_avg=measures_avg, selected=selected))
}


#' @export
opnmfR_ranksel_perm <- function(X, rs, W0=NULL, use.rcpp=TRUE, nperm=1, plots=TRUE, seed=NA, rtrue=NA, ...) {
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
    if(!is.na(rtrue)) abline(v=rtrue, lty=2, col="gray")
    legend("bottomleft", legend = c("Orig.","Perm."), pch=16:17)
    title(main="Permutation", sub=paste("Selected rank = ", min(selr)))
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
  # the following conditions might not be met when doing ooser
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

