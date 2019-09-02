#' X =  d x N matrix
#' U = d x k matrix
r1pca.Cr <- function(X,U,cc, xxt){
  N <- ncol(X)
  d <- nrow(X)

  if(nrow(U) != d)
    stop("U Matrix rows does not match that of X")

  message("Computing w")
  resid <- X - U %*% t(U) %*% X
  w <- sqrt(apply(resid^2,2,sum))
  w <- ifelse(w <= cc, 1, cc/w)

  message("Computing Cr")
  ans <- matrix(0,nrow=d,ncol=d)
  for(i in 1:N){
    ans <- ans + w[i] * xxt[[i]]
  }

  ans
}

r1pca.residues <- function(X, U){
  N <- ncol(X)
  d <- nrow(X)

  if(nrow(U) != d)
    stop("U Matrix rows does not match that of X")

  s <- sapply(1:N,function(i)
    sqrt(t(X[,i])%*%X[,i] - t(X[,i])%*%U%*%t(U)%*%X[,i]))

  s
}

r1pca <- function(X,k,niter = 100){
  N <- ncol(X)
  U <- irlba::irlba(X, nv = k)$u
  s <- r1pca.residues(X,U)
  cc <- median(s)

  message("Calculating xxt")
  xxt <- list()
  xxt <- lapply(1:N, function(i) X[,i] %*% t(X[,i]))

  for(i in 1:niter){
    message(paste0("Error: ",sum(r1pca.residues(X,U))))
    Cr <- r1pca.Cr(X,U,cc, xxt)
    U <- Cr%*%U

    message("Performing gram schmidt")
    U <- pracma::gramSchmidt(U)$Q
  }
  list(u = U, 
       d = t(U)%*%Cr%*%U,
       v = t(U)%*%X)
}

