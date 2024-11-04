library(cglasso)

### setting ####
n <- 8
L <- 3
nstars.size <- 3
nstars <- 1
p <- nstars.size * nstars * (L + 1)
q <- 3
n.ThtY <- p * ( p - 1 ) / 2
n.b <- p*q


ntp <- tmp <- 30

a <- 1
B.min <- a * 1.4 # imput parameter for rPar()
B.max <- a * 2.1 # imput parameter for rPar()

b <- 1
tht.min <- b * 2.4  # imput parameter for rPar()
tht.max <- b * 2.5  # imput parameter for rPar()

g.ebic <- 0.5

Omg <- diag(q)*1
scal.terms <- c(1,0.2) # ensure decreasing traces of the covariance matrices of the scores over L
scal.terms[1]*c(1:L)^{-scal.terms[2]}

nsim <- 60
niter <- 3

rho.min.ratio <- 0.05
rho.max.ratio <- 1
nrho <- 20
perc.rho.seq <- seq(from = rho.max.ratio , to= rho.min.ratio , length.out = nrho)
lmbd.min.ratio <- 0.05
lmbd.max.ratio <- 1
nlambda <- 20
perc.lmbd.seq <- seq(from= lmbd.max.ratio , to= lmbd.min.ratio , length.out = nlambda)

B.supp <- matrix(0, nrow = q, ncol = p,
                 dimnames = list(paste0("X", seq_len(q)),
                                 paste0("Y", seq_len(p))))
set.seed(1)
for(k in seq_len(p))
  B.supp[, k] <- sample(c(0,1,1))

Par <- rPar(p = p, q = q, L = L, B.supp = B.supp, B.min = B.min, 
            B.max = B.max, nstars.size = nstars.size, nstars = nstars,  
            tht.min = tht.min, tht.max = tht.max, Omg = Omg, scal.terms = scal.terms, union = F )
B.true <- vector( mode = "list", length = L)
for( l in seq_len(L)){
  sig.x <- Par$Sgm[1:q,1:q,l]
  B.true[[l]] <-  solve(sig.x) %*% Par$Sgm[1:q, -c(1:q),l]
  print(diag(Par$Sgm[-c(1:q),-c(1:q),l] - t(B.true[[l]]) %*% sig.x %*% B.true[[l]]))
}

Tht.supp <- matrix(0, p, p)
for (l in seq_len(L)){
  Tht.supp <- Tht.supp + (Par$Tht[ - seq_len(q), - seq_len(q), l] != 0)
}

Tht.supp[upper.tri(Tht.supp, TRUE)] <- 0
E.tht <- which(Tht.supp != 0, arr.ind = TRUE)

###### Tht and Sgm [Y X]
Sgm.YX <- Tht.YX <- array(NA , dim = dim(Par$Tht))
colnames(Sgm.YX) <- colnames(Tht.YX) <- rownames(Tht.YX) <- rownames(Sgm.YX) <- c(paste0("Y", seq_len(p)), paste0("X", seq_len(q)))

for (l in seq_len(L)){
  Sgm.YX[ 1:p , 1:p , l] <- Par$Sgm[-c(1:q),-c(1:q),l]
  Sgm.YX[ -c(1:p) , -c(1:p) , l] <- Par$Sgm[1:q,1:q,l]
  Sgm.YX[ 1:p , -c(1:p) , l] <- Par$Sgm[-c(1:q),1:q,l]
  Sgm.YX[ -c(1:p) , 1:p , l] <- Par$Sgm[1:q,-c(1:q),l]
  Tht.YX[ , ,l ] <- solve(Sgm.YX[ , , l])
}
Sigma_trueY <- array(NA, dim = c(p,p,L))
for(l in seq_len(L)) Sigma_trueY[,,l] <- Par$Sgm[-(1:q),-(1:q) , l]

weights <- rep(n, L) / (n * L)

b <- c(round(range(Par$B), 2), round(range(Par$Tht[-c(1:q), -c(1:q), ]), 2 ))
#####################

lmbd_estim <- array(NA, dim = c(nsim, 6, niter) , dimnames = list( 
  "sim" = paste0("sim_", seq(from = 1, to = nsim)), c("Oracle KL", "Oracle Acc YB","Oracle Acc B","KLCV", "AIC", "eBIC") , 
  "iteration" = paste0( "iter_", seq(1, niter))))
rho_estim <- array(NA, dim = c(nsim, 5, niter) , dimnames = list( 
  "sim" = paste0("sim_", seq(from = 1, to = nsim)), c("Oracle KL", "Oracle Acc","KLCV", "AIC", "eBIC") , 
  "iteration" = paste0( "iter_", seq(1, niter))))

true_KL_loss <- array(NA, dim = c(nsim, 4, niter, 2), 
                      dimnames = list( "sim" = paste0("sim_", seq(from = 1, to = nsim)),
                                       "KL_loss" = c("Oracle" ,"KLCV", "AIC", "eBIC"),
                                       "iter" = paste0("iter_", seq(from = 1, to = niter)),
                                       c("ThetaYX", "ThetaY") ))
Accuracy <- array(NA, dim = c(nsim, 4, niter, 3), 
                  dimnames = list( "sim" = paste0("sim_", seq(from = 1, to = nsim)),
                                   "Acc" = c("Oracle","KLCV", "AIC", "eBIC"),
                                   "iter" = paste0("iter_", seq(from = 1, to = niter)),
                                   c("ThetaYX", "ThetaY", "B") ))
AUC <- array(NA, dim = c(nsim, 4, niter, 3), 
             dimnames = list( "sim" = paste0("sim_", seq(from = 1, to = nsim)),
                              "AUC" = c("Oracle","KLCV", "AIC", "eBIC"),
                              "iter" = paste0("iter_", seq(from = 1, to = niter)),
                              c("B", "ThetaY.lmbd.rho", "ThetaY.rho") ))
oracle.last.iter <- matrix(NA, nsim, 2, dimnames = list("sim" = paste0("sim_", seq(from = 1, to = nsim)), c("KL.Y", "KL.YX")) )

tht.stimYX.klcv <- tht.stimYX.AIC <- tht.stimYX.eBIC <- vector( mode = "list", length = L)
for (l in seq_len(L)) {
  tht.stimYX.klcv[[l]] <- tht.stimYX.AIC[[l]] <- tht.stimYX.eBIC[[l]] <- array( NA, dim = c( p+q, p+q, 1) )
}

pb <- txtProgressBar(min = 0L, max = nsim, style = 3L)
jj <- 0L
set.seed(1)
for(h in seq_len(nsim)) {
  jj <- jj + 1L
  setTxtProgressBar(pb, jj)
  df <- rcfggm(n = n, p = p, q = q, ntp = ntp, Sgm = Par$Sgm)$Xi.z
  
  x <- y <-  df1 <- Sgm <- data.list <- eps.Y <- B.stim <- out <- vector(mode = "list", length = L)
  l <- 1
  for( l in seq_len(L)) {
    x[[l]] <- df[,1:q,l]
    y[[l]] <- df[,-c(1:q),l]
    Sgm[[l]] <- Par$Sgm[-c(1:q), -c(1:q), l]
    
    data <- datacggm( Y = y[[l]], X = x[[l]] )
    out[[l]] <- cglasso(.~. , data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
    
    eps.Y[[l]] <- matrix(NA, n, q)
    B.stim[[l]] <- matrix(NA, p, q)# array(NA, dim = c(q, p, nrho) )
    eps.Y[[l]] <- y[[l]]- out[[l]]$mu[,,2,1]
    B.stim[[l]] <- t( out[[l]]$B[-1, , 2, ])
  }
  data.list <- datajcggm(Y = y, X = x)
  
  ### rho
  rho_max <- jcglasso(data = datajcggm(Y = eps.Y), nrho = 1, alpha1 = 0)$rho
  rho1 <- perc.rho.seq * rho_max
  out.tht <-  jcglasso(data = datajcggm(Y = eps.Y), rho = rho1, alpha1 = 0, nu = 0, alpha3 = 0)

  tht.stim <- sig.stim <- vector(mode = "list", length = L)
  l <- 1
  for(l in seq_len(L))
  {
    tht.stim[[l]] <- out.tht$Tht[ out.tht$InfoStructure$id_Y, out.tht$InfoStructure$id_Y, l, 1 , ] 
    sig.stim[[l]] <- out.tht$Sgm[out.tht$InfoStructure$id_Y, out.tht$InfoStructure$id_Y, l, 1, ]
  }

  id.min.klcv <- which.min( KLCV.fixed.X(tht.stim = tht.stim, eps.Y = eps.Y, L ) )
  id.min.aic <- which.min(AIC(out.tht)$value_gof)
  id.min.ebic <- which.min(BIC(out.tht, g = .5, type = "FD" )$value_gof)
  
  rho_estim[jj, "KLCV", 1] <- out.tht$rho[ id.min.klcv ]
  rho_estim[jj, "AIC", 1] <- out.tht$rho[ id.min.aic ]
  rho_estim[jj, "eBIC", 1] <- out.tht$rho[ id.min.ebic ]
  
  KL_loss <- KL.tht1(Theta_stim = tht.stim , Theta_true = Par$Tht[-c(1:q), -c(1:q), ], L = L,
                     npar = nrho, X = x, B_stim = B.stim, B_true = Par$B )
  
  true_KL_loss[ jj , "Oracle", 1, "ThetaY"] <- min(KL_loss)
  rho_estim[jj, "Oracle KL", 1] <- out.tht$rho[which.min(KL_loss)]
  true_KL_loss[ jj , "KLCV", 1, "ThetaY"] <- KL_loss[ id.min.klcv ]
  true_KL_loss[ jj , "AIC", 1, "ThetaY"] <- KL_loss[ id.min.aic ]
  true_KL_loss[ jj , "eBIC", 1, "ThetaY"] <- KL_loss[ id.min.ebic ]
  
  thetah <- array(0, dim = c(p,p,nrho))
  for (l in seq_len(L)) thetah <- thetah + tht.stim[[l]]
  accuracy <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho)$accuracy
  rho_estim[jj, "Oracle Acc", 1] <-  out.tht$rho[which.max(accuracy)]
  Accuracy[ jj , "Oracle", 1, "ThetaY"] <- max(accuracy)
  Accuracy[ jj , "KLCV", 1, "ThetaY"] <- accuracy[ id.min.klcv ]
  Accuracy[ jj , "AIC", 1, "ThetaY"] <- accuracy[ id.min.aic ]
  Accuracy[ jj , "eBIC", 1, "ThetaY"] <- accuracy[ id.min.ebic ]
  
  AUC[ jj, "Oracle", 1, "ThetaY.rho"] <- AUC[ jj, "KLCV", 1, "ThetaY.rho"] <- AUC[ jj, "AIC", 1, "ThetaY.rho"] <- AUC[ jj, "eBIC", 1, "ThetaY.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
  
  rm(id.min.klcv, id.min.aic, id.min.ebic, accuracy, KL_loss, tht.stim, sig.stim, out.tht, rho_max, rho1)
  
  print( rbind( "KL Y " = round(true_KL_loss[ jj , , 1, "ThetaY"],3), "acc Y " = round( Accuracy[ jj , , 1, "ThetaY"], 3)), "AUC" =   AUC[ jj, 1, 1, "ThetaY.rho"] )
  
  ##### iterations start
  iter <- 1
  for ( iter in seq_len(niter-1)) {
    #### oracle KL
    lambda.oracle.kl <- lambda.iteration( data = data.list , rho = rho_estim[jj, "Oracle KL" , iter], L = L, p = p, 
                                          perc.lmbd.seq = perc.lmbd.seq, weights = weights)
    
    KLTRUE <- KL.tht_B(ThetaXY_stim = lambda.oracle.kl$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q) 

    lmbd_estim[jj, "Oracle KL", iter ] <- lambda.oracle.kl$lambda[ which.min( KLTRUE )]
    true_KL_loss[ jj , "Oracle", iter, "ThetaYX"] <- min( KLTRUE )
    
    rm( lambda.oracle.kl , KLTRUE )
    #### oracle ACC 
    lambda.oracle.acc <- lambda.iteration( data = data.list , rho = rho_estim[jj, "Oracle Acc" , iter], L = L,  p = p,
                                           perc.lmbd.seq = perc.lmbd.seq, weights = weights)
    thetah <- array(0, dim = c(p,p,nlambda))
    bh <- array(0, dim = c(q,p,nlambda))
    for (l in seq_len(L)) {
      thetah <- thetah + lambda.oracle.acc$tht.stimY[[l]]
      bh <- bh + round(lambda.oracle.acc$b.stim[[l]],10)
    }
    AccuracyTHTY <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                        accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                        ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                        accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                       ( n.ThtY + n.b ))
    Accuracy[ jj ,"Oracle", iter, "ThetaYX"] <- max( AccuracyTHTY )      
    lmbd_estim[jj, "Oracle Acc YB", iter ] <- lambda.oracle.acc$lambda[which.min( AccuracyTHTY )]
    
    AccuracyB <- ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )
    Accuracy[jj , "Oracle", iter, "B"] <- max( AccuracyB )
    lmbd_estim[jj, "Oracle Acc B", iter ] <- lambda.oracle.acc$lambda[which.min( AccuracyB )]
    
    AUC[ jj, "Oracle", iter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
    AUC[ jj, "Oracle", iter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    rm(lambda.oracle.acc , thetah, bh, AccuracyTHTY, AccuracyB )
    
    ### lambda KLCV
    lambda.klcv <- lambda.iteration( data = data.list , rho = rho_estim[jj, "KLCV" , iter], L = L, p = p,
                                     perc.lmbd.seq = perc.lmbd.seq, weights = weights)
    id.min.klcv <- which.min( KLCV.gauss2(tht.stimYX = lambda.klcv$tht.stimYX, b.stim = lambda.klcv$b.stim, data = lambda.klcv$data, L = L, whole.tht = F ))
    lmbd_estim[jj, "KLCV", iter ] <- lambda.klcv$lambda[ id.min.klcv ]
    true_KL_loss[ jj , "KLCV", iter, "ThetaYX"] <-  KL.tht_B(ThetaXY_stim = lam <- lambda.klcv$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[ id.min.klcv ]

    thetah <- array(0, dim = c(p,p,nlambda))
    bh <- array(0, dim = c(q,p,nlambda))
    for (l in seq_len(L)) {
      thetah <- thetah + lambda.klcv$tht.stimY[[l]]
      bh <- bh + round(lambda.klcv$b.stim[[l]],10)
    }
    Accuracy[ jj , "KLCV", iter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                                                   accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                                   ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                                                   accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                  ( n.ThtY + n.b ))[id.min.klcv]
    
    Accuracy[jj , "KLCV", iter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                             accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.klcv]
    
    AUC[ jj, "KLCV", iter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
    AUC[ jj, "KLCV", iter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    rm(lambda.klcv , id.min.klcv , thetah, bh )
    ### lambda AIC
    lambda.aic <- lambda.iteration( data = data.list , rho = rho_estim[jj, "AIC", iter], L = L, p = p,
                                    perc.lmbd.seq = perc.lmbd.seq, weights = weights)
    id.min.aic <- which.min( AIC.jcglasso(lambda.aic$out.jcglasso)$value_gof )
    lmbd_estim[jj, "AIC", iter ] <- lambda.aic$lambda[ id.min.aic ] 
    true_KL_loss[ jj , "AIC", iter, "ThetaYX"] <- KL.tht_B(ThetaXY_stim = lambda.aic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[ id.min.aic ]

    thetah <- array(0, dim = c(p,p,nlambda))
    bh <- array(0, dim = c(q,p,nlambda))
    for (l in seq_len(L)) {
      thetah <- thetah + lambda.aic$tht.stimY[[l]]
      bh <- bh + round(lambda.aic$b.stim[[l]],10)
    }
    Accuracy[ jj , "AIC", iter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TP +
                                                  ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                 ( n.ThtY + n.b ))[id.min.aic]
    
    Accuracy[jj , "AIC", iter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                            accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.aic]
    
    AUC[ jj, "AIC", iter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
    AUC[ jj, "AIC", iter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    rm(lambda.aic , id.min.aic, thetah, bh )
    
    ### lambda eBIC
    lambda.ebic <- lambda.iteration( data = data.list , rho = rho_estim[jj,"eBIC", iter], L = L,  p = p, 
                                     perc.lmbd.seq = perc.lmbd.seq, weights = weights)
    id.min.ebic <- which.min(BIC.jcglasso(lambda.ebic$out.jcglasso, g = g.ebic, type = "FD")$value_gof )
    lmbd_estim[jj, "eBIC", iter ] <- lambda.ebic$lambda[id.min.ebic]
    true_KL_loss[ jj , "eBIC", iter, "ThetaYX"] <- KL.tht_B(ThetaXY_stim = lambda.ebic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[ id.min.ebic]
    
    thetah <- array(0, dim = c(p,p,nlambda))
    bh <- array(0, dim = c(q,p,nlambda))
    for (l in seq_len(L)) {
      thetah <- thetah + lambda.ebic$tht.stimY[[l]]
      bh <- bh + round(lambda.ebic$b.stim[[l]],10)
    }
    Accuracy[ jj , "eBIC", iter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                                                   accuracy.B.2( bh = bh, B.true = B.supp)$TP +
                                                   ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                                                   accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                  ( n.ThtY + n.b ) )[id.min.ebic]
    
    Accuracy[jj , "eBIC", iter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                             accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.ebic]
    
    AUC[ jj, "eBIC", iter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
    AUC[ jj, "eBIC", iter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    rm(lambda.ebic , id.min.ebic, thetah, bh )
    
    print( rbind( "KL YX " = round(true_KL_loss[ jj , , iter, "ThetaYX"],3), 
                  "acc YX " = round( Accuracy[ jj , , iter, "ThetaYX"], 3),
                  "AUC B " = round( AUC[ jj , , iter, "B"], 3),
                  "AUC Y " = round( AUC[ jj , , iter, "ThetaY.lmbd.rho"], 3)))
    
    #
    #
    
    ### rho oracle KL
    rho.oracle.kl <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[jj, "Oracle KL", iter], L = L, perc.rho.seq = perc.rho.seq)
    KLTRUE <- KL.tht_B(ThetaXY_stim = rho.oracle.kl$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q)

    rho_estim[jj, "Oracle KL", iter + 1] <-  rho.oracle.kl$rho[which.min(  KLTRUE )]
    true_KL_loss[ jj , "Oracle", iter + 1 , "ThetaY"] <- min(KLTRUE)
    
    rm(rho.oracle.kl, KLTRUE )
    ### rho oracle acc
    rho.oracle.acc <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[jj, "Oracle Acc YB", iter], 
                                     L = L, perc.rho.seq = perc.rho.seq)
    thetah <- array(0, dim = c(p,p,nrho))
    for (l in seq_len(L)) thetah <- thetah + rho.oracle.acc$tht.stim[[l]][1:p, 1:p,]
    accuracy.rho <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho)$accuracy
    Accuracy[ jj , "Oracle", iter+1 , "ThetaY"] <- max(accuracy.rho)
    rho_estim[jj, "Oracle Acc", iter + 1] <- rho.oracle.acc$rho[which.max(accuracy.rho)]
    
    AUC[jj, "Oracle", iter +1 , "ThetaY.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    rm(rho.oracle.acc,   accuracy.rho , thetah )
    
    #### rho KLCV
    rho.klcv <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[jj, "KLCV", iter], L = L, perc.rho.seq = perc.rho.seq)
    id.min.klcv <- which.min( KLCV.gauss2( tht.stimYX = rho.klcv$tht.stim, b.stim = rho.klcv$b.stim, data = rho.klcv$data, L = L, whole.tht = F ) )

    rho_estim[jj, "KLCV", iter + 1] <- rho.klcv$rho[id.min.klcv ]
    true_KL_loss[ jj , "KLCV", iter+1 , "ThetaY"] <- KL.tht_B(ThetaXY_stim = rho.klcv$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q)[ id.min.klcv]

    thetah <- array(0, dim = c(p,p,nrho))
    for (l in seq_len(L)) thetah <- thetah + rho.klcv$tht.stim[[l]][1:p, 1:p,]
    Accuracy[ jj , "KLCV", iter+1 , "ThetaY"] <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho)$accuracy[ id.min.klcv ]
    
    AUC[ jj, "KLCV", iter + 1, "ThetaY.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    if( iter == niter - 1 ) oracle.last.iter[jj, "KL.Y"] <- min(KL.tht_B(ThetaXY_stim = rho.klcv$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q))

    rm(rho.klcv , id.min.klcv, thetah )
    
    #### rho AIC
    rho.aic <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[jj, "AIC", iter], L = L, perc.rho.seq = perc.rho.seq)
    id.min.aic <- which.min(AIC(rho.aic$out)$value_gof)
    rho_estim[jj, "AIC", iter + 1] <- rho.aic$rho[id.min.aic ]
    true_KL_loss[ jj , "AIC", iter +1, "ThetaY"] <- KL.tht_B(ThetaXY_stim = rho.aic$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q)[id.min.aic ]

    thetah <- array(0, dim = c(p,p,nrho))
    for (l in seq_len(L)) thetah <- thetah + rho.aic$tht.stim[[l]][1:p, 1:p,]
    Accuracy[ jj , "AIC", iter+1 , "ThetaY"] <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho)$accuracy[ id.min.aic ]
    
    AUC[ jj, "AIC", iter + 1, "ThetaY.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    a <- min(KL.tht_B(ThetaXY_stim = rho.aic$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q))
    if( iter == niter - 1 &  a < oracle.last.iter[jj, "KL.Y"] ) oracle.last.iter[jj, "KL.Y"] <- a
    
    rm(rho.aic , id.min.aic, thetah, a )
    
    #### rho eBIC
    rho.ebic <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[jj,  "eBIC", iter], L = L, perc.rho.seq = perc.rho.seq)
    id.min.ebic <- which.min(BIC(rho.ebic$out, g = .5, type = "FD")$value_gof)
    rho_estim[jj, "eBIC", iter + 1] <- rho.ebic$rho[id.min.ebic]
    true_KL_loss[ jj , "eBIC", iter + 1, "ThetaY"] <- KL.tht_B(ThetaXY_stim = rho.ebic$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q)[ id.min.ebic]
    thetah <- array(0, dim = c(p,p,nrho))
    for (l in seq_len(L)) thetah <- thetah + rho.ebic$tht.stim[[l]][1:p,1:p,]
    Accuracy[ jj , "eBIC", iter+1 , "ThetaY"] <- ROCTheta(thetah = thetah, E = E.tht, nrho = nrho)$accuracy[ id.min.ebic ]
    
    AUC[ jj, "eBIC", iter + 1, "ThetaY.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
    
    b <- min(KL.tht_B(ThetaXY_stim = rho.ebic$tht.stim, ThetaXY_true = Par$Tht, L = L, npar= nrho , p = p, q = q))
    if( iter == niter - 1 &  b < oracle.last.iter[jj, "KL.Y"] ) oracle.last.iter[jj, "KL.Y"] <- b
    
    rm(rho.ebic , id.min.ebic, thetah, b )
    
    print( rbind( "KL Y " = round(true_KL_loss[ jj , , iter, "ThetaY"],3), "acc Y " = round( Accuracy[ jj , , iter, "ThetaY"], 3)))
    
    iter <- iter + 1
  }
  
  ################## end for{ niterations }
  # final estimate for [ chi gamma ]
  
  #### oracle KL
  lambda.oracle.kl <- lambda.iteration( data = data.list , rho = rho_estim[jj, "Oracle KL" , niter], L = L,  p = p,
                                        perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  KLTRUE <- KL.tht_B(ThetaXY_stim = lambda.oracle.kl$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)

  lmbd_estim[jj, "Oracle KL", niter ] <- lambda.oracle.kl$lambda[ which.min(KLTRUE) ]
  true_KL_loss[ jj , "Oracle", niter, "ThetaYX"] <- min(KLTRUE)
  rm( lambda.oracle.kl , KLTRUE )
  #### oracle ACC 
  lambda.oracle.acc <- lambda.iteration( data = data.list , rho = rho_estim[jj, "Oracle Acc" , niter], L = L, p = p,
                                         perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  thetah <- array(0, dim = c(p,p,nlambda))
  bh <- array(0, dim = c(q,p,nlambda))
  for (l in seq_len(L)) {
    thetah <- thetah + lambda.oracle.acc$tht.stimY[[l]]
    bh <- bh + round(lambda.oracle.acc$b.stim[[l]],10)
  }
  AccuracyTHTY <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                      accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                      ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                      accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                     ( n.ThtY + n.b ))
  Accuracy[ jj ,"Oracle", niter, "ThetaYX"] <- max( AccuracyTHTY )      
  lmbd_estim[jj, "Oracle Acc YB", niter ] <- lambda.oracle.acc$lambda[which.min( AccuracyTHTY )]
  
  AccuracyB <- ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )
  Accuracy[jj , "Oracle", niter, "B"] <- max( AccuracyB )
  lmbd_estim[jj, "Oracle Acc B", niter ] <- lambda.oracle.acc$lambda[which.min( AccuracyB )]
  
  
  AUC[ jj, "Oracle", niter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
  AUC[ jj, "Oracle", niter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
  
  
  rm(lambda.oracle.acc , thetah, bh, AccuracyTHTY, AccuracyB )
  
  ### lambda KLCV
  lambda.klcv <- lambda.iteration( data = data.list , rho = rho_estim[jj, "KLCV", niter], L = L, p = p,
                                   perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  id.min.klcv <- which.min( KLCV.gauss2(tht.stimYX = lambda.klcv$tht.stimYX, b.stim = lambda.klcv$b.stim, data = lambda.klcv$data, L = L, whole.tht = F ))
  lmbd_estim[jj, "KLCV", niter] <- lambda.klcv$lambda[id.min.klcv]
  true_KL_loss[ jj , "KLCV", niter, "ThetaYX"] <- KL.tht_B(ThetaXY_stim = lambda.klcv$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[id.min.klcv]
  thetah <- array(0, dim = c(p,p,nlambda))
  bh <- array(0, dim = c(q,p,nlambda))
  for (l in seq_len(L)) {
    thetah <- thetah + lambda.klcv$tht.stimY[[l]]
    bh <- bh + round(lambda.klcv$b.stim[[l]],10)
  }
  Accuracy[ jj , "KLCV", niter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP + 
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TP +
                                                  ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN + 
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                 ( n.ThtY + n.b  ) )[id.min.klcv]
  Accuracy[jj , "KLCV", niter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                            accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.klcv]
  
  
  
  AUC[ jj, "KLCV", niter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
  AUC[ jj, "KLCV", niter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
  
  oracle.last.iter[jj, "KL.YX"] <- min(KL.tht_B(ThetaXY_stim = lambda.klcv$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q))

  rm(lambda.klcv , id.min.klcv, thetah, bh)
  
  ### lambda AIC
  lambda.aic <- lambda.iteration( data = data.list , rho = rho_estim[jj, "AIC", niter], L = L, p = p,
                                  perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  id.min.aic <- which.min( AIC.jcglasso(lambda.aic$out)$value_gof )
  lmbd_estim[jj, "AIC", niter] <- lambda.aic$lambda[id.min.aic]
  true_KL_loss[ jj , "AIC", niter, "ThetaYX"] <- KL.tht_B(ThetaXY_stim = lambda.aic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[ id.min.aic]

  thetah <- array(0, dim = c(p,p,nlambda))
  bh <- array(0, dim = c(q,p,nlambda))
  for (l in seq_len(L)) {
    thetah <- thetah + lambda.aic$tht.stimY[[l]]
    bh <- bh + round(lambda.aic$b.stim[[l]],10)
  }
  Accuracy[ jj , "AIC", niter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP +
                                                 accuracy.B.2( bh = bh, B.true = B.supp)$TP +
                                                 ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN +
                                                 accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                ( n.ThtY + n.b ) )[id.min.aic]
  
  Accuracy[jj , "AIC", niter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                           accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.aic]
  
  
  AUC[ jj, "AIC", niter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
  AUC[ jj, "AIC", niter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
  
  a <- min( KL.tht_B(ThetaXY_stim = lambda.aic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q))
  if( a < oracle.last.iter[jj, "KL.YX"] )   oracle.last.iter[jj, "KL.YX"] <- a
  
  rm(lambda.aic , id.min.aic, thetah, bh, a)
  
  ### lambda eBIC
  lambda.ebic <- lambda.iteration( data = data.list , rho = rho_estim[jj,"eBIC", niter], L = L, p = p,
                                   perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  id.min.ebic <- which.min(BIC.jcglasso(lambda.ebic$out, g = g.ebic, type = "FD")$value_gof )
  lmbd_estim[jj, "eBIC", niter] <- lambda.ebic$lambda[id.min.ebic]
  
  true_KL_loss[ jj , "eBIC", niter, "ThetaYX"] <- KL.tht_B(ThetaXY_stim = lambda.ebic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q)[id.min.ebic]
  thetah <- array(0, dim = c(p,p,nlambda))
  bh <- array(0, dim = c(q,p,nlambda))
  for (l in seq_len(L)) {
    thetah <- thetah + lambda.ebic$tht.stimY[[l]]
    bh <- bh + round(lambda.ebic$b.stim[[l]],10)
  }
  Accuracy[ jj , "eBIC", niter, "ThetaYX"] <- ((ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TP +
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TP +
                                                  ROCTheta(thetah = thetah, E = E.tht, nrho = nlambda)$TN +
                                                  accuracy.B.2( bh = bh, B.true = B.supp)$TN) /
                                                 ( n.ThtY + n.b ) )[id.min.ebic]
  
  Accuracy[jj , "eBIC", niter, "B"] <-  ((accuracy.B.2( bh = bh, B.true = B.supp)$TP + 
                                            accuracy.B.2( bh = bh, B.true = B.supp)$TN) / n.b )[id.min.ebic]
  
  
  AUC[ jj, "eBIC", niter, "B"] <- ROCB(Bh = bh, Bt = B.supp, nlambda = nlambda)$AUC
  AUC[ jj, "eBIC", niter, "ThetaY.lmbd.rho"] <- ROCTheta(thetah = thetah, E = E.tht , nrho = nrho)$ROC
  
  b <- min( KL.tht_B(ThetaXY_stim = lambda.ebic$tht.stimYX, ThetaXY_true = Par$Tht, L = L, npar= nlambda , p = p, q = q))
  if( b < oracle.last.iter[jj, "KL.YX"] )   oracle.last.iter[jj, "KL.YX"] <- b
  
  
  rm(lambda.ebic , id.min.ebic, thetah, bh, b)
  # iteration completed for rho and lambda
  
  if( niter == 1){
    if(jj > 1) print( rbind( "median KL YX " = round(apply( true_KL_loss[1:jj,,,"ThetaYX"], 2, median),3), "median acc YX " = round( apply( Accuracy[1:jj,,,"ThetaYX"], 2, median), 3)))
    else print( rbind( "median KL YX " = round(true_KL_loss[jj,,,"ThetaYX"], 3), "median acc YX " = round( Accuracy[jj,,,"ThetaYX"], 3)))
  }
  
  if( niter != 1){
    if(jj > 1) print( rbind( "median KL YX " = round(apply( true_KL_loss[1:jj,,,"ThetaYX"], c(2,3) , median),3), "median acc YX " = round( apply( Accuracy[1:jj,,,"ThetaYX"], c(2,3), median), 3)))
    else print( rbind( "median KL YX " = round(true_KL_loss[jj,,,"ThetaYX"], 3), "median acc YX " = round( Accuracy[jj,,,"ThetaYX"], 3)))
  }
  
}


# ORACLE
median(oracle.last.iter[, "KL.Y"])
median(oracle.last.iter[, "KL.YX"])
median(Accuracy[ , "Oracle" ,niter,"ThetaY"])
median(Accuracy[ , "Oracle" ,niter,"B"])
round(apply(AUC[,"Oracle", niter, ], 2, median), 2)


# KLCV AIC EBIC
# YX
round(apply( true_KL_loss[1:26,,niter,"ThetaYX"], 2, median),2) # apply( true_KL_loss[,,niter,"ThetaYX"], 2, mean)
round(apply( Accuracy[ 1:26,,niter,"ThetaYX"], 2, median),2) #; apply( Accuracy[,,niter,"ThetaYX"], 2, mean)
# Y
round(apply( true_KL_loss[ ,,niter,"ThetaY"], 2, median),2) # apply( true_KL_loss[,,niter,"ThetaY"], 2, mean)
round(apply( Accuracy[ ,,niter,"ThetaY"], 2, median),2) #; apply( Accuracy[,,niter,"ThetaY"], 2, mean)
# B
round(apply( Accuracy[ ,,niter,"B"], 2, median),2) #; apply( Accuracy[,,niter,"ThetaB"], 2, mean)

round(apply( Accuracy[,-1,niter,"ThetaYX"], 2, median),2); round(apply( Accuracy[,-1,niter,"ThetaYX"], 2, mean),2)
round(apply( Accuracy[,-1,niter,"ThetaY"], 2, median),2); round(apply( Accuracy[,-1,niter,"ThetaY"], 2, mean),2)
apply( Accuracy[,-1,niter,"B"], 2, median); apply( Accuracy[,-1,niter,"B"], 2, mean)








