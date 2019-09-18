NBrotate <- function(X, y, V0 = NULL, ndim = NULL, method = 'Gaussian', hmult = 1){
  n <- nrow(X)
  d <- ncol(X)
  yy <- y
  u <- unique(y)
  for(k in 1:length(u)) yy[which(y==u[k])] = k
  if(is.null(ndim)){
    if(is.null(V0)) ndim <- max(yy)-1
    else ndim <- ncol(V0)
  }
  else if(ndim > d) stop('cannot compute projection with more dimensions than the data')
  PI <- sapply(1:length(u), function(k) sum(yy==k))/n
  MU <- t(sapply(1:max(yy), function(k) colMeans(X[which(yy==k),])))
  if(is.null(V0)){
    M <- sapply(1:max(yy), function(k) yy==k)%*%MU
    Sw <- cov(X-M)
    Sb <- cov(X)-Sw
    V0 <- eigen(solve(Sw + diag(d)*1e-10)%*%Sb + cov(X)*1e-10)$vectors[,1:ndim]
  }
  S <- lapply(1:length(u), function(k) cov(X[which(yy==k),]))
  if(method=='Gaussian'){
    V <- optim(V0, fnb_G, dfnb_G, X, yy, MU, S, PI, n, method = 'BFGS', control = list(fnscale = -1))$par
    V <- sweep(V, 2, apply(V, 2, norm_vec), '/')
    sol <- list(V=V, X = X, y = yy, MU=MU, S=S, PI=PI, nc = max(yy), y_labels = u, method = 'Gaussian')
    class(sol) <- "NBrotate"
    sol
  }
  else if(method=='kernel'){
    mn <- colMeans(X)
    X <- sweep(X, 2, mn, '-')
    h <- hmult*sapply(1:length(u), function(k) mean(eigen(S[[k]])$values[1:ndim])^.5/(n*PI[k])^(1/5))
    V <- optim(V0, f_NB_2, df_NB_2, X, yy, PI, h, method = 'BFGS', control = list(fnscale = -1))$par
    V <- sweep(V, 2, apply(V, 2, norm_vec), '/')
    sol <- list(V=V, X = X, y=yy, PI=PI, h=h, nc = max(yy), mu = mn, y_labels = u, method = 'kernel')
    class(sol) <- "NBrotate"
    sol
  }
  else stop('"method" must be one of "Gaussian" and "kernel"')
}




predict.NBrotate <- function(sol, Xte){
  if(sol$method=='Gaussian'){
    posterior <- matrix(0,nrow(Xte),sol$nc)
    colnames(posterior) <- sol$y_labels
    mu <- sol$MU%*%sol$V
    s <- lapply(sol$S, function(s) diag(t(sol$V)%*%s%*%sol$V))
    for(k in 1:sol$nc) posterior[,k] <- exp(-distmat(Xte%*%t(t(sol$V)/sqrt(s[[k]])),mu[k,]/sqrt(s[[k]]))^2/2)/prod(s[[k]]^.5)*sol$PI[k]
    posterior <- posterior[,order(sol$y_labels)]
    list(class = sort(sol$y_labels)[apply(posterior, 1, which.max)], posterior = posterior/rowSums(posterior))
  }
  else if(sol$method=='kernel'){
    nte <- nrow(Xte)
    dd <- ncol(sol$V)
    denc <- matrix(0,nte,sol$nc)
    p <- sol$X%*%sol$V
    pte <- sweep(Xte, 2, sol$mu, '-')%*%sol$V
    ote <- apply(pte, 2, order)
    ixs <- lapply(1:sol$nc, function(k) which(sol$y==k))
    oc <- lapply(ixs, function(l) apply(p[l,], 2, order))
    for(c in 1:sol$nc){
      den_c <- numeric(nte) + 1
      for(dm in 1:dd){
        den_c[ote[,dm]] <- den_c[ote[,dm]]*den_fast(p[ixs[[c]][oc[[c]][,dm]],dm],pte[ote[,dm],dm],sum(sol$y==c),nte,sol$h[c])
      }
      denc[,c] <- sol$PI[c]*den_c
    }
    posterior <- denc/rowSums(denc)
    posterior <- posterior[,order(sol$y_labels)]
    colnames(posterior) <- sort(sol$y_labels)
    list(class = sort(sol$y_labels)[apply(posterior, 1, which.max)], posterior = posterior)
  }
}
