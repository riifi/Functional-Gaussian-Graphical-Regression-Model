source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/datajcggm.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/jcglasso.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/jcggm.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/gof.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/to_graph.R")

download.file("https://github.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/raw/main/02%20-%20Simulations/SimulationStudy2/CGLASSO/cglasso_2.0.5.tar.gz", destfile = "cglasso_2.0.5.tar.gz")
#install.packages("cglasso_2.0.5.tar.gz", repos = NULL, type = "source")

rPar <- function(p, q, L, B.supp, B.min, B.max, nstars.size, nstars, tht.min = 1 / 3, tht.max = 2 / 3, Omg = diag(q), scal.terms = c(1,1.8), union = F ) {
  #--- Setting Output
  K <- p + q
  B <- array(0, dim = c(q, p, L),
             dimnames = list(paste0("X", seq_len(q)),
                             paste0("Y", seq_len(p)),
                             paste0("L", seq_len(L))))
  Sgm <- Tht <- array(0, dim = c(K, K, L), 
                      dimnames = list(c(paste0("X", seq_len(q)), paste0("Y", seq_len(p))),
                                      c(paste0("X", seq_len(q)), paste0("Y", seq_len(p))),
                                      paste0("L", seq_len(L))))
  #--- sampling matrix B
  for(k in seq_len(p)) {
    id <- B.supp[, k] == 1
    nid <- sum(id)
    B[id, k, ] <- runif(L * nid, B.min, B.max) * (rbinom(L * nid, 1, prob = 0.5) * 2 - 1)
  }
  #-- sampling covariance and precision matrix
  for (l in seq_len(L)) {
    if(union)
      #El.id <- c("c", paste0("L", l))
      #El <- E[E.id %in% El.id, ]
      Tht[-seq_len(q), -seq_len(q), l] <- rStar.l(nstars.size = nstars.size, 
                                                  nstars = nstars, L = L, l = l, 
                                                  tht.min = tht.min, 
                                                  tht.max = tht.max)
    if(!union) Tht[-seq_len(q), -seq_len(q), l] <- rStar.l1(nstars.size = nstars.size, 
                                                            nstars = nstars, L = L,
                                                            tht.min = tht.min, 
                                                            tht.max = tht.max)
    
    Tht[seq_len(q), -seq_len(q), l] <- - B[, , l] %*% Tht[-seq_len(q), -seq_len(q), l]
    Tht[-seq_len(q), seq_len(q), l] <- t(Tht[seq_len(q), -seq_len(q), l])
    Tht[seq_len(q), seq_len(q), l] <- Omg + B[, , l] %*% Tht[-seq_len(q), -seq_len(q), l] %*% t(B[, , l])
    Tht[, , l] <- Tht[, , l] / (scal.terms[1] * l^( - scal.terms[2]))
    Sgm[, , l] <- solve(Tht[, , l])
  }
  return(list(B = B, Tht = Tht, Sgm = Sgm))
}

rcfggm <- function(n, p, q, ntp, Sgm){
  # ------------
  # Arguments
  #
  # n: sample size
  # p: number of response curves
  # q: number of covariate curves
  # ntp: number of time points
  # Sgm: array di dimensione (p + q) x (p + q) x L
  # sgm2.eps: error variance
  # ------------
  L <- dim(Sgm)[3L]
  tp <- seq(from = 0, to = 1, length = ntp)
  # OUTPUT
  Z <- array(0, dim = c(ntp, n, dim(Sgm)[1L]),
             dimnames = list(paste0("tp.", seq_len(ntp)),
                             NULL,
                             paste0("Z", seq_len(dim(Sgm)[1L]))))
  # Starting Sim
  Xi.z <- array(0, dim = c(n, dim(Sgm)[1L], L), 
                dimnames = list(1:n, 
                                paste0("Z", seq_len(dim(Sgm)[1L])),
                                paste0("l", seq_len(L))))
  for (l in seq_len(L)) {
    Xi.z[, , l] <- MASS::mvrnorm(n = n, 
                                 mu = rep(0, dim(Sgm)[1L]),
                                 Sigma = Sgm[, , l])
  }
  return(list(Z = Z, tp = tp, Xi.z = Xi.z))
}

KLCV.fixed.X <- function(tht.stim, eps.Y, L ){
  npar <- dim(tht.stim[[1]])[3]
  N <- dim( eps.Y[[l]])[1]
  KLCV <- rep(0, npar)
  S <- vector(mode = "list", length = L)
  for(r in seq_len(npar)){
    for(l in seq_len(L)){
      S <- t( eps.Y[[l]])%*% eps.Y[[l]]/N
      KLCV[r] <- KLCV[r] - (1/(2*L) )* (log(det(tht.stim[[l]][,,r])) - sum(diag(tht.stim[[l]][,,r] %*% S)))
      
      # BIAS
      I <- tht.stim[[l]][,,r] != 0
      bias <- rep(NA, nrho)
      for(n in seq_len(N)){
        Sn <-  eps.Y[[l]][n,] %*% t( eps.Y[[l]][n,])
        KLCV[r] <- KLCV[r] + (1/(L*N*(N*L-1))) * sum((solve(tht.stim[[l]][,,r]) - Sn) * I * (tht.stim[[l]][,,r]%*%((S-Sn)* I) %*% tht.stim[[l]][,,r] ))
      }}}
  KLCV <- KLCV
  return(KLCV) 
}

KL.tht1 <- function(Theta_stim, Theta_true, L, npar, X, B_stim, B_true){ # Sigma's dimensions pxpxL ----> Par$Sgm[ p+q , p+q , L] 
  # out$Tht$class1[ p , p , 1, nrho]
  Sigma_true <- array(NA, dim = dim(Theta_true))
  for( l in seq_len(L)) Sigma_true[ , , l] <- solve(Theta_true[, , l])
  p <- dim(Sigma_true)[1]
  N <- dim(X[[1]])[1]
  KL <- vector(mode = "double", length = npar)
  for (k in seq_len(npar)) {
    for (l in seq_len(L)) {
      A <- Sigma_true[ , , l ] %*% Theta_stim[[l]][ , , k ] #[ , , ,k ]
      KL[k] <- KL[k]+ 1/2 * ( sum(diag( Theta_stim[[l]][ , , k ] %*% (( t(B_true[,,l]) - B.stim[[l]]) %*% t(X[[l]]) %*% X[[l]] %*% t( t( B_true[,,l] ) - B.stim[[l]])) ))/N + sum(diag( A )) - log( det( A ) ) - p)
    }
  }
  return(KL)
}

ROCTheta <- function(thetah, E, nrho){
  p <- dim(thetah)[1]
  Tht.sup <- matrix(0, p, p)
  Tht.sup[E] <- 1
  U <- lower.tri(Tht.sup, diag = FALSE)
  thetat_v <- Tht.sup[U]
  A <- which(abs(thetat_v) > 0)
  thetah_m <- apply(thetah, 3, function(M) M[U])
  P.rho <- apply(thetah, 3, function(M) ( sum(M!=0) -p )/2 )
  
  TPR <- if(length(A) == 1) apply(t(abs(thetah_m[A, ]) > 0), 2, mean) else apply(abs(thetah_m[A, ]) > 0, 2, mean)
  TP <- if(length(A) == 1) apply(t(abs(thetah_m[A, ]) > 0), 2, sum) else apply(abs(thetah_m[A, ]) > 0, 2,sum)
  FPR <- vector(mode = "numeric", length = nrho)
  
  for(i in 1:nrho){
    id <- which(abs(thetah_m[, i]) > 0)
    FP <- sum(!(id %in% A))
    TN <- sum(which(thetah_m[, i] == 0) %in% which(thetat_v == 0))
    FPR[i] <- FP/(FP+TN)
  }
  
  TPR.roc <- c(0, TPR, 1)
  FPR.roc <- c(0, FPR, 1)
  
  dTPR <- c(diff(TPR.roc), 0)
  dFPR <- c(diff(FPR.roc), 0)
  ROC <- sum(na.omit(TPR.roc * dFPR)) + sum(dTPR * dFPR) / 2
  
  class.corr <- apply(thetah_m, 2, function(M) sum(as.numeric(M !=0) == thetat_v))/length(thetat_v)
  
  out <- list(TPR = TPR, FPR = FPR, ROC = ROC, TN = TN, TP = TP,P.rho =  P.rho, accuracy = class.corr)
  out
}

lambda.iteration <- function( data , rho, L, perc.lmbd.seq, p, weights ){
  
  tmp <- jcglasso(data = data, rho = rho, lambda = 1E10, nu = 0, 
                  alpha1 = 0, alpha2 = 0, alpha3 = 0)
  
  lambda_max <- max(sqrt(rowSums(sapply(seq_len(L), function(k) {
    R <- tmp$R[,,k,1,1] / n
    X <- tmp$Zipt[,tmp$InfoStructure$id_X,k,1,1]
    Tht <- coef(tmp, "Theta", class.id = k, rho.id = 1, lambda.id = 1)
    (weights[k] * (t(X) %*% R %*% Tht))^2
  }))))
  
  lambda1 <- perc.lmbd.seq * lambda_max
  
  out <- jcglasso(data = data , rho = rho, lambda = lambda1, nu = 0, 
                  alpha1 = 0, alpha2 = 0, alpha3 = 0)
  
  b.stim <- tht.stimYX <- tht.stimY <- sig.stim <- tht.stimY <- vector(mode = "list", length = L)
  for(l in seq_len(L))
  {
    b.stim[[l]] <- out$B[-1, , l, , ]
    tht.stimYX[[l]] <- out$Tht[,, l, , ]
    tht.stimY[[l]] <- tht.stimYX[[l]][ 1:p , 1:p ,]
  }
  
  out <- list( lambda = lambda1 , out.jcglasso = out, tht.stimYX = tht.stimYX, tht.stimY = tht.stimY,
               b.stim = b.stim, data = out$Z )
  return( out )
}

rho.iteration <- function( y, x, n, p, q, lambda, L, perc.rho.seq ){
  
  tmp <- jcglasso(data = data.list, rho = 1E10, lambda = lambda,, nu = 0, 
                  alpha1 = 0, alpha2 = 0, alpha3 = 0)
  
  rho_max <- max(sqrt(rowSums(sapply(seq_len(L), function(k) {
    y <- tmp$Zipt[,tmp$InfoStructure$id_Y,k,1,1]
    mu <- fitted(tmp, class.id = k, lambda.id = 1)
    R <- y - mu
    (weights[k] * crossprod(R) / n)[outer(1:p, 1:p, "<")]^2
  }))))
  
  rho1 <- perc.rho.seq * rho_max
  out.tht <-  jcglasso(data = data.list, rho = rho1, lambda = lambda , nu = 0, 
                       alpha1 = 0, alpha2 = 0, alpha3 = 0)
  
  tht.stim <- sig.stim <- b <- vector(mode = "list", length = L)
  for(l in seq_len(L))
  {
    tht.stim[[l]] <- out.tht$Tht[ , , l, 1 , ] 
    sig.stim[[l]] <- out.tht$Sgm[ , , l, 1, ]
    b[[l]] <- out.tht$B[ -1, ,l ,1 , ]
  }
  out <- list( rho = out.tht$rho , out = out.tht,  
               tht.stim = tht.stim , b.stim = b , data = out.tht$Z)
  return(out)
}

rho.iteration.first <- function( y, x, n, p, q, lambda, L, perc.rho.seq ){
  for( l in seq_len(L)) {
    data <- datacggm( Y = y[[l]], X = x[[l]] )
    out[[l]] <- cglasso(.~. , data = data, lambda =  lambda, nrho = 1)
    eps.Y[[l]] <- matrix(NA, n, q)
    B.stim[[l]] <- matrix(NA, p, q)# array(NA, dim = c(q, p, nrho) )
    eps.Y[[l]] <- datajcggm( Y =( y[[l]] - out[[l]]$mu[,,1,1] ))
    B.stim[[l]] <- t( out[[l]]$B[-1, , 1, ])
  }
  rho_max <- jcglasso(data = eps.Y, nrho = 1, alpha = 0)$rho
  rho1 <- perc.rho.seq * rho_max
  out.tht <-  jcglasso(data = eps.Y, rho = rho1, alpha = 0)
  
  tht.stim <- sig.stim <- vector(mode = "list", length = L)
  Sigma_true <- array(NA, dim = c(p,p,L))
  for(l in seq_len(L))
  {
    tht.stim[[l]] <- out.tht$Tht[ , , l, 1 , ] 
    sig.stim[[l]] <- out.tht$Sgm[ , , l, 1, ]
    Sigma_true[,,l] <- Par$Sgm[-(1:q),-(1:q) , l]
  }
  out <- list( rho = out.tht$rho , out.jcglasso = out.tht,  out.cglasso = out,  
               tht.stim = tht.stim, B.stim = B.stim , data = eps.Y)
  return(out)
}

KL.tht_B <- function(ThetaXY_stim, ThetaXY_true, L, npar, p, q){ 
  SigmaY_true <- array(NA, dim = c(p,p,L) )
  SigmaX_true <- array(NA, dim = c(q,q,L) )
  for( l in seq_len(L)){
    SigmaX_true[ , , l] <- solve(ThetaXY_true[1:q,1:q , l])
    SigmaY_true[ , , l] <- solve(ThetaXY_true[-(1:q),-(1:q), l])
  }
  KL <- vector(mode = "double", length = npar)
  for (k in seq_len(npar)) {
    for (l in seq_len(L)) {
      XA <-  SigmaX_true[ , , l ] %*% ThetaXY_stim[[l]][-(1:p),-(1:p), k ] 
      YA <- SigmaY_true[ , , l ] %*% ThetaXY_stim[[l]][1:p,1:p, k ] #[ , , ,k ]
      KL[k] <- KL[k]+ 1/2 * (sum(diag( XA )) - log( det( XA ) ) - q + sum(diag( YA )) - log( det( YA ) ) - p)
    }
  }
  return(KL)
}

ROCB <- function(Bh, Bt, nlambda){
  Bt_v <- c(Bt)
  A <- which(abs(Bt_v) > 0)
  Bh_m <- apply(Bh, 3, function(M) c(M))
  
  TPR <- if(length(A) == 1) apply(t(abs(Bh_m[A, ]) > 0), 2, mean) else apply(abs(Bh_m[A, ]) > 0, 2, mean)
  FPR <- vector(mode = "numeric", length = nlambda)
  
  for(i in 1:nlambda){
    id <- which(abs(Bh_m[, i]) > 0)
    FP <- sum(!(id %in% A))
    TN <- sum(which(Bh_m[, i] == 0) %in% which(Bt_v == 0))
    FPR[i] <- FP/(FP+TN)
  }
  TPR.roc <- c(0, TPR, 1)
  FPR.roc <- c(0, FPR, 1)
  
  dTPR <- c(diff(TPR.roc), 0)
  dFPR <- c(diff(FPR.roc), 0)
  AUC <- sum(na.omit(TPR.roc * dFPR)) + sum(dTPR * dFPR) / 2
  out <- list(TPR = TPR, FPR = FPR, AUC = AUC)
  out
}

KLCV.gauss2 <- function(tht.stimYX, b.stim, data, L, whole.tht ){
  npar <- dim(tht.stimYX[[1]])[3]
  N <- dim(data[[l]]$Y)[1]
  KLCV <- rep(0, npar)
  S.b <- vector(mode = "list", length = L)
  #  H <- length(b.stim)     H and L are the same
  for(r in seq_len(npar))
  {
    for(l in seq_len(L))
    {
      if( !whole.tht){
        S.b.y <- t(data[[l]]$Y - data[[l]]$X %*% b.stim[[l]][,,r])%*%(data[[l]]$Y - data[[l]]$X %*% b.stim[[l]][,,r])/N
        S.x <- (t(data[[l]]$X) %*% data[[l]]$X) / N
        if( is.matrix(tht.stimYX[[l]][-(1:p),-(1:p),r] ) ) 
          KLCV[r] <- KLCV[r] - ( 1/(2*L) ) * ( log(det(tht.stimYX[[l]][-(1:p),-(1:p),r])) - sum(diag(tht.stimYX[[l]][-(1:p),-(1:p),r] %*% S.x)) + log(det(tht.stimYX[[l]][1:p,1:p,r])) - sum(diag(tht.stimYX[[l]][1:p,1:p,r] %*% S.b.y)))
        else 
          KLCV[r] <- KLCV[r] - ( 1/(2*L) ) * ( - tht.stimYX[[l]][-(1:p),-(1:p),r] %*% S.x + log(det(tht.stimYX[[l]][1:p,1:p,r])) - sum(diag(tht.stimYX[[l]][1:p,1:p,r] %*% S.b.y)))
      }
      
      if( whole.tht ){
        eps.y <- data[[l]]$Y - data[[l]]$X %*% b.stim[[l]][,,r]
        S <- (t(cbind(eps.y, data[[l]]$X)) %*% cbind(eps.y, data[[l]]$X)) / N
        dimnames(S)[[1]] <- dimnames(S)[[2]] <- dimnames(tht.stimYX[[1]])$response
        KLCV[r] <- KLCV[r] - ( 1/(2*L) ) * ( log(det(tht.stimYX[[l]][,,r])) - sum(diag(tht.stimYX[[l]][,,r] %*% S)))
        
      }

      # BIAS
      if( !whole.tht){
        I.x <- tht.stimYX[[l]][-(1:p),-(1:p),r] != 0
        I.y <- tht.stimYX[[l]][1:p,1:p,r] != 0
        bias <- rep(NA, nrho)
        for(n in seq_len(N)){
          Sn.b.y <- t(data[[l]]$Y[n,] - data[[l]]$X[n,] %*% b.stim[[l]][,,r])%*%(data[[l]]$Y[n,] - data[[l]]$X[n,] %*% b.stim[[l]][,,r])
          Sn.x <- data[[l]]$X[n,]%*%t(data[[l]]$X[n,])
          
          a <- sum( (solve(tht.stimYX[[l]][-(1:p),-(1:p),r]) - Sn.x) * I.x * (tht.stimYX[[l]][-(1:p),-(1:p),r]%*%((S.x - Sn.x) * I.x) %*% tht.stimYX[[l]][-(1:p),-(1:p),r] ))
          b <- sum( (solve(tht.stimYX[[l]][1:p,1:p,r]) - Sn.b.y) * I.y * (tht.stimYX[[l]][1:p,1:p,r]%*%((S.b.y-Sn.b.y)* I.y) %*% tht.stimYX[[l]][1:p,1:p,r] ))
          KLCV[r] <- KLCV[r] + 1/(L*N*(N*L-1)) * (a + b)
        }
      }
      if(whole.tht){
        I <- tht.stimYX[[l]][,,r] != 0
        bias <- rep(NA, nrho)
        for(n in seq_len(N)){
          eps.y.n <- data[[l]]$Y[n,] - data[[l]]$X[n,] %*% b.stim[[l]][,,r]
          Sn <- c(eps.y.n, data[[l]]$X[n,])%*%t(c(eps.y.n, data[[l]]$X[n,]))
          dimnames(Sn)[[1]] <- dimnames(Sn)[[2]] <- dimnames(tht.stimYX[[1]])$response
          a <- sum( (solve(tht.stimYX[[l]][,,r]) - Sn) * I * ( tht.stimYX[[l]][, , r]%*%((S - Sn) * I) %*% tht.stimYX[[l]][,,r] ))
          KLCV[r] <- KLCV[r] + 1/(L*N*(N*L-1)) * a
        }
      }
    }
  }
  KLCV
  return(KLCV) 
}

accuracy.B.2 <- function(bh, B.true){
  p <- dim(B.true)[2]
  q <- dim(B.true)[1]
  npar <- dim(bh)[3]
  Pos <- which(B.true > 0)
  Neg <- which(B.true == 0)
  Pos.h <- apply(bh, 3, function(M) which( M != 0 ))
  if( sum( bh != 0) == 0) Pos.h <- as.list( rep( 0, npar) )
  Neg.h <- apply(bh, 3, function(M) which( M == 0 ))
  if( sum( bh == 0) == 0) Neg.h <- as.list( rep( 0, npar) )
  
  TP <- TN <- vector(mode = "integer", length = npar)
  for( j in seq_len(npar)){
    TP[j] <- sum(Pos.h[[j]] %in% Pos)
    TN[j] <- sum(Neg.h[[j]] %in% Neg)
  }
  
  accuracy2 <- ( TP/length(Pos) + TN/length(Neg) ) - 1
  
  out <- list(accuracy2 = accuracy2, TP = TP, TN = TN, P = length(Pos), N = length(Neg))
  out
}

rStar.l1 <- function(nstars.size, nstars = 1, L, tht.min = 0.4, tht.max = 0.5) {
  # Arguments:
  # nstars.size: number of vertices of each star
  # nstars: number of stars in each block
  # L: total number of layers (the number of block is L + 1 as the first block 
  #     is common to all layers)
  # l: 
  # tht.min, tht.max: the same as in rStar.R
  require("pracma")
  Tht.l <- diag(nstars.size * nstars * L)
  blk <- rep(seq_len(L), each = nstars.size * nstars)
  for (l in seq_len(L)) {
    Tht.l[blk == l, blk == l] <- rStars(d = nstars.size, nstars = nstars, 
                                        tht.min = tht.min, tht.max = tht.max) 
    
  }
  
  Tht.c <- rStars(d = nstars.size, nstars = nstars, tht.min = tht.min, tht.max = tht.max)
  return(blkdiag(Tht.c, Tht.l))
}

rStars <- function(d, nstars = 1, tht.min = 0.4, tht.max = 0.5) {
  # Arguments
  # d: number of outer vertices 
  # nstars: number of stars
  # tht.min, tht.max:
  require("pracma")
  out <- vector(mode = "list", length = nstars)
  for(l in seq_len(nstars)) {
    Tht <- diag(d)
    Tht.yx <- runif(d - 1, tht.min, tht.max) * sample(c(-1, 1), d - 1, TRUE)
    K <- outer(Tht.yx, Tht.yx)
    Tht.xx <- diag((1 + sqrt(1 + 4 * diag(K))) / 2)
    Tht[1, 1] <- 1 + drop(Tht.yx %*% solve(Tht.xx) %*% Tht.yx)
    Tht[1, -1] <- Tht[-1, 1] <- Tht.yx
    Tht[-1, -1] <- Tht.xx
    out[[l]] <- Tht
  }
  return(do.call(blkdiag, out))
}

