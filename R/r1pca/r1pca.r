#' X =  d x N matrix
#' U = d x k matrix
r1pca.Cr <- function(X,U,cc){
  N <- ncol(X)
  d <- nrow(X)

  if(nrow(U) != d)
    stop("U Matrix rows does not match that of X")

  resid <- X - U %*% t(U) %*% X
  w <- sqrt(apply(resid^2,2,sum))
  if(any(w == 0))
    stop("One of the ws = 0, error in projection")
  w <- ifelse(w<=cc, 1, cc/w)

  X <- t(t(X)*sqrt(w))
  X %*% t(X)
}

r1pca.residues <- function(X, U){
  N <- ncol(X)
  d <- nrow(X)

  if(nrow(U) != d)
    stop("U Matrix rows does not match that of X")

  ans <- sapply(1:N,function(i)
    sqrt(t(X[,i])%*%X[,i] - t(X[,i])%*%U%*%t(U)%*%X[,i]))

  # Debug
  if(any(is.nan(ans)))
    print(ans)

  ans
}

r1pca <- function(X,k,verbose  = F){
  N <- ncol(X)
  U <- irlba::irlba(X, nv = k)$u
  s <- r1pca.residues(X,U)
  cc <- median(s)

  prev.err <- Inf
  cur.err <- sum(r1pca.residues(X,U))
  while(prev.err > cur.err){
    prev.err <- cur.err

    if(verbose)
      message(paste0("Error = ",cur.err))
    Cr <- r1pca.Cr(X,U,cc)
    U <- Cr%*%U
    U <- qr.Q(qr(U))

    cur.err <- sum(r1pca.residues(X,U))
    if(is.nan(cur.err) | is.na(cur.err))
      cur.err <- Inf
  }
  list(u = U, 
       d = t(U)%*%Cr%*%U,
       v = t(U)%*%X)
}

