norm_vec <- function(v) sqrt(sum(v^2))


f_NB_2 <- function(V, X, y, mypi, h){
  V <- matrix(V,nrow=ncol(X))
  n <- nrow(X)
  d <- ncol(X)
  dd <- ncol(V)
  nc <- length(h)
  #den <- array(0,dim = c(n,dd,nc)) might not be needed for function evaluation, only derivative
  den <- matrix(1, n, nc)
  p <- t(t(X%*%V)/apply(V,2,norm_vec))
  o <- apply(p, 2, order)
  ns <- sapply(1:max(y), function(i) sum(y==i))
  for(dm in 1:dd){
    den[o[,dm],] <- den[o[,dm],]*den_fast_all(p[o[,dm],dm], y[o[,dm]]-1, n, ns, nc, max(ns), h)
  }
  den[which(is.na(den))] <- 1e-300
  den[which(den<=0)] <- 1e-300
  T1 <- sum(sapply(1:nc, function(i) sum(log(mypi[i]*den[which(y==i),i]))))
  T2 <- sum(log(rowSums(sweep(den, 2, mypi, '*'))))
  T1-T2
}



df_NB_2 <- function(V, X, y, mypi, h){
  V <- matrix(V,nrow=ncol(X))
  n <- nrow(X)
  d <- ncol(X)
  dd <- ncol(V)
  nc <- length(h)
  #den <- array(0,dim = c(n,dd,nc)) might not be needed for function evaluation, only derivative
  den <- array(0, dim = c(n, dd, nc))
  p <- t(t(X%*%V)/apply(V,2,norm_vec))
  o <- apply(p, 2, order)
  ns <- sapply(1:max(y), function(i) sum(y==i))
  for(dm in 1:dd){
    den[o[,dm],dm,] <- den_fast_all(p[o[,dm],dm], y[o[,dm]]-1, n, ns, nc, max(ns), h)
  }
  den[which(is.na(den))] <- 1e-300
  den[which(den<=0)] <- 1e-300
  denc <- apply(den, c(1,3), prod)
  pr <- rowSums(sweep(denc, 2, mypi, '*'))
  dp <- p*0
  for(dm in 1:dd){
    dencc <- denc/den[,dm,]
    dp[o[,dm],dm] <- dp[o[,dm],dm] + dksum_all_in_class(p[o[,dm],dm], matrix(1, n, nc), y[o[,dm]]-1, n, ns, nc, max(ns), h)*(1/n/mypi/h^2)[y[o[,dm]]]/den[,dm,][cbind(o[,dm],y[o[,dm]])]
    dp[o[,dm],dm] <- dp[o[,dm],dm] + dksum_all_in_class(p[o[,dm],dm], 1/den[o[,dm],dm,], y[o[,dm]]-1, n, ns, nc, max(ns), h)*(1/n/mypi/h^2)[y[o[,dm]]]
    dp[o[,dm],dm] <- dp[o[,dm],dm] - rowSums(sweep(dksum_cross(p[o[,dm],dm], matrix(1, n, nc), y[o[,dm]]-1, n, ns, nc, max(ns), h)/pr[o[,dm]], 2, n*h^2, '/')*dencc[o[,dm],])
    dp[o[,dm],dm] <- dp[o[,dm],dm] - dksum_all(p[o[,dm],dm], sweep(dencc/pr, 2, mypi/ns/h^2, '*')[o[,dm],], y[o[,dm]]-1, n, nc, h)
  }
  dV <- V*0
  for(i in 1:dd) dV[,i] <- dp[,i]%*%(X/norm_vec(V[,i])-p[,i]%*%t(V[,i])/norm_vec(V[,i])^2)
  c(dV)
}



NB_classify <- function(V, X, Xte, y, mypi, h){
  V <- matrix(V,nrow=ncol(X))
  n <- nrow(X)
  nte <- nrow(Xte)
  d <- ncol(X)
  dd <- ncol(V)
  nc <- length(h)
  denc <- matrix(0,nte,nc)
  p <- t(t(X%*%V)/apply(V,2,norm_vec))
  pte <- t(t(Xte%*%V)/apply(V,2,norm_vec))
  ote <- apply(pte, 2, order)
  oc <- lapply(ixs, function(l) apply(p[l,], 2, order))
  for(c in 1:nc){
    den_c <- numeric(nte) + 1
    for(dm in 1:dd){
      den_c[ote[,dm]] <- den_c[ote[,dm]]*den_fast(p[ixs[[c]][oc[[c]][,dm]],dm],pte[ote[,dm],dm],sum(y==c),nte,h[c])
    }
    denc[,c] <- pi[c]*den_c
  }
  posterior <- denc/rowSums(denc)
  list(class = apply(posterior, 1, which.max), posterior = posterior)
}
