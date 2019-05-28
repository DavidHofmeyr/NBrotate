### function evaluation for projection matrix V
## V = projection matrix
## X = data matrix
## Y = vector of class labels
## MU = matrix of class means
## S = list of class covariance matrices
## PI = vector of class priors
## N = total number of data (i.e., nrow(X))

fnb_G <- function(V, X, Y, MU, S, PI, N){
  ### normalise V and turn into a matrix in case it is given as a vector
  V <- matrix(V, nrow = ncol(MU))
  nv <- apply(V, 2, norm_vec)
  v <- t(t(V)/nv)

  ### compute statistics on projection
  mu <- MU%*%v
  s <- lapply(S, function(s) diag(t(v)%*%s%*%v))


  nc <- length(PI)

  ### compute projected data (in each class. This could be removed if we always use the actual covariances in S,
  ### or scaled versions of it. It might be useful to try other alternatives, though, such as adding a ridge
  ### term to decrease the flexibility)
  x <- lapply(1:nc, function(k) X[which(Y==k),]%*%v)


  ### The first term in the log-likelihood essentially only depends on the variance terms:
  T1 <- sum(N*PI*log(PI)-sapply(1:nc, function(k) N*PI[k]*sum(log(s[[k]]))/2+(N*PI[k]-1)*sum(diag(cov(x[[k]]))/s[[k]])))

  ### compute the (scaled) density of all points in all classes
  p <- matrix(0,N,nc)
  for(k in 1:nc) p[,k] <- exp(-distmat(X%*%t(t(v)/sqrt(s[[k]])),mu[k,]/sqrt(s[[k]]))^2/2)/prod(s[[k]]^.5)*PI[k]


  T1 - sum(log(rowSums(p)))
}

dfnb_G <- function(V, X, Y, MU, S, PI, N){
  ### The first steps in the gradient are the same as in the function evaluation

  V <- matrix(V, nrow = ncol(MU))
  nc <- length(PI)
  nv <- apply(V, 2, norm_vec)
  v <- t(t(V)/nv)
  mu <- MU%*%v
  sv <- lapply(S, function(s) s%*%v)
  s <- lapply(S, function(s) diag(t(v)%*%s%*%v))
  p <- matrix(0,nrow(X),nc)
  for(k in 1:nc) p[,k] <- exp(-distmat(X%*%t(t(v)/sqrt(s[[k]])),mu[k,]/sqrt(s[[k]]))^2/2)/prod(s[[k]]^.5)*PI[k]
  ps <- rowSums(p)

  ### The derivative of the first term in the log-lkelihood also only depends on the derivatives of the variance terms
  ssv <- lapply(1:nc, function(k) t(t(sv[[k]])/s[[k]]))
  dv1 <- v*0
  for(k in 1:nc) dv1 <- dv1 - ssv[[k]]*N*PI[k]

  ### Next compute the derivative of the second term in the log-likelihood
  dv2 <- v*0
  for(k in 1:nc){
    Xk <- sweep(X, 2, MU[k,], '-')
    dv2 <- dv2 + sweep(sv[[k]], 2, colSums((p[,k]/ps)*sweep(sweep(Xk%*%v, 2, sqrt(s[[k]]), '/')^2-1, 2, s[[k]], '/')), '*')
    dv2 <- dv2 - t((p[,k]/ps)*Xk)%*%Xk%*%t(t(v)/s[[k]])
  }
  dv <- dv1-dv2
  for(k in 1:ncol(V)) dv[,k] <- dv[,k]/nv[k] - t(dv[,k])%*%v[,k]%*%t(v[,k])/nv[k]
  dv
}
