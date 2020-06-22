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
opnmfR <- function(X,r,W0=NULL,max.iter=50000,tol=1e-5,memsave=TRUE,eps=1e-16,...) {
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
    
    W[W < eps] <- eps
    W <- W / norm(W,"2") # spectral norm
    
    setTxtProgressBar(pbar, iter)
    
    diffW <- norm(Wold-W, "F") / norm(Wold, "F")
    if(diffW < tol) break
  }
  close(pbar)
  
  
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

#' @export
opnmfR_ranksel_perm <- function(X, rs, W0=NULL, use.rcpp=TRUE, plots=TRUE, seed=NA, ...) {
  stopifnot(ncol(X)>=max(rs))
  start_time <- Sys.time()
  
  if(is.na(seed)) seed <- sample(1:10^6, 1)
  set.seed(seed)
  
  Xperm <- apply(X,2,sample) # permute the rows in each column
  nn <- list()
  mse <- list()
  for(i in 1:length(rs)) {
    if(use.rcpp) {
      nn[[i]] <- opnmfRcpp(X, rs[i], W0=W0, ...)
      nnp <- opnmfRcpp(Xperm, rs[i], W0=W0, ...)
    } else {
      nn[[i]] <- opnmfR(X, rs[i], W0=W0, ...)
      nnp <- opnmfR(Xperm, rs[i], W0=W0, ...)
    }
    
    cat("orig", "iter:", nn[[i]]$iter, "time:", nn[[i]]$time, "diffW:", nn[[i]]$diffW, "\n")
    cat("perm", "iter:", nnp$iter, "time:", nnp$time,"diffW:", nnp$diffW, "\n")
    
    mse[[i]] <- list()
    mse[[i]]$orig <- opnmfR_mse(X, nn[[i]]$W, nn[[i]]$H)
    mse[[i]]$perm <- opnmfR_mse(Xperm, nnp$W, nnp$H)
  }
  
  names(mse) <- rs
  names(nn) <- rs
  
  mseorig <- sapply(mse, function(xx) xx$orig)
  mseperm <- sapply(mse, function(xx) xx$perm)
  mseorig <- (mseorig-min(mseorig)) / (max(mseorig)-min(mseorig))
  mseperm <- (mseperm-min(mseperm)) / (max(mseperm)-min(mseperm))
  
  sel <- which(diff(mseorig) > diff(mseperm))
  if(length(sel)>1) sel <- sel
  selr <- rs[sel]
  
  if(plots) {
    maxi <- max(c(mseorig, mseperm))
    mini <- min(c(mseorig, mseperm))
    
    plot(NA, xlim=c(min(rs),max(rs)), ylim=c(mini,maxi), xlab="Rank", ylab="MSE (scaled)")
    points(rs, mseorig, type='b', pch=16)
    points(rs, mseperm, type='b', pch=17, lty=2)
    points(selr, mseorig[sel], cex=2, pch=1, lwd=2, col="red")
    legend("topright", legend = c("Orig.","Perm."), pch=16:17)
    title(main="Permutation based rank selection", sub=paste("Selected rank = ", selr))
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

