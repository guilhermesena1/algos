###% Auxliary functions
.pos <- function(x){ as.numeric(x>=0) * x }
.neg <- function(x){ - as.numeric(x<0) * x }
.norm <- function(x){ sqrt(drop(crossprod(x))) }

# Initial guess of W and H via nonnegative SVD
NMF.InitialGuess <- function(A, k){
	#check the input matrix
	if( any(A<0) ) stop('The input matrix contains negative elements !')
	
	#size of input matrix
	size = dim(A);
	m <- size[1]; n<- size[2]
	
	#the matrices of the factorization
	W = matrix(0, m, k);
	H = matrix(0, k, n);
	
	#1st SVD --> partial SVD rank-k to the input matrix A.	
	s = svd(A, k, k, LINPACK=LINPACK);
	U <- s$u; S <- s$d; V <- s$v
	
	#choose the first singular triplet to be nonnegative
	W[,1] = sqrt(S[1]) * abs(U[,1]);         
	H[1,] = sqrt(S[1]) * abs(t(V[,1])); 
					
	# second SVD for the other factors (see table 1 in Boutsidis' paper)
	for( i in seq(2,k) ){
		uu = U[,i]; vv = V[,i];
		uup = .pos(uu); uun = .neg(uu) ;
		vvp = .pos(vv); vvn = .neg(vv);
		n_uup = .norm(uup);
		n_vvp = .norm(vvp) ;
		n_uun = .norm(uun) ;
		n_vvn = .norm(vvn) ;
		termp = sum(n_uup %*% n_vvp); 
    termn = sum(n_uun %*% n_vvn);
		if (termp >= termn){
			W[,i] = sqrt(S[i] * termp) * uup / n_uup; 
			H[i,] = sqrt(S[i] * termp) * vvp / n_vvp;
		}else{		
			W[,i] = sqrt(S[i] * termn) * uun / n_uun; 
			H[i,] = sqrt(S[i] * termn) * vvn / n_vvn;
		}
	}
	
	#------------------------------------------------------------
	#actually these numbers are zeros
	W[W<0.0000000001] <- 0;
	H[H<0.0000000001] <- 0;

  W <- as(W, "sparseMatrix")
  H <- as(H, "sparseMatrix")

  return (list(W = W, H = H))
}

# Runs NMF with or without masking (masking estimates CV error)
nmf <- function(A, k = 10, err = 1e-1, cross = F, mask = NULL){
  A2 <- NULL
  cve <- NULL

  # Zeros some values of A to impute later
  if(cross){
    A2 <- A
    # If mask not given choose 10% of A elements randomly
    if(is.null(mask)){
      mask <- sample(length(A), 0.1*length(A))
    }
    
    A[which(mask == T)] <- 0
  }

  init <- NMF.InitialGuess(A, k)
  W <- init$W
  H <- init$H
  rm(init)

  prev <- Inf
  cur <- sum((A - W%*%H)^2)

  while((!is.nan(cur)) & (prev - cur > err)){
    Wt <- t(W)
    Ht <- t(H)

    # Update H
    tmph1 <- Wt %*% A;
    tmph2 <- Wt %*% W %*% H
    H <- H * tmph1/tmph2

    # Update W
    tmpw1 <- A %*% Ht
    tmpw2 <- W %*% H %*% Ht
    W <- W * tmpw1/tmpw2

    prev <- cur
    cur <- sum((A - W%*%H)^2)
  }

  #CV : find imputation MSE
  if(cross){
    B <- W%*%H
    cve <- sum((A2[mask] - B[mask])^2)
  }
  list(W = W, H = H, cve = cve)
}

# Tries different values of k and cross-validates
cv.nmf <- function(A, ks = 1:50, fast.cv = F){
  tenfold <- cut(1:length(A), 10)
  tenfold <- sample(tenfold)
  vals <- levels(tenfold)
  cves <- sapply(
    foreach(k = ks) %dopar% {
      ans <- 0
      if(fast.cv) {
        vals <- vals[1:2]
      }
      for(val in vals) {
        mask <- (tenfold == val)
        res <- nmf(A, k, cross = T, mask = mask)
        ans <- ans + res$cve
      }
      print(paste0("CVE[",k,"] = ",ans))
      ans
    }, sum)

  cves <- cves/length(A)
  names(cves) <- paste0("k_", ks)
  return (cves)
}
