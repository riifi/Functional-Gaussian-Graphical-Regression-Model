##########
## TO RUN IF YOU DOWNLOADED THE FILE "IAGOS smooth data.RData" FROM
## https://github.com/riifi/Functional-Gaussian-Graphical-Regression-Model/tree/main
## "IAGOS smooth data.RData" CONTAINS SMOOTHED DATA ( y_n and x_n as exposed in the section "Evaluation of the scores")
##
## If you prefer to download the data from https://iagos.aeris-data.fr/download/ to get
## the smoothed data from the original dataset, the "IAGOS data import.R" file contains the code to import 
## the nc4-format data and generate the scores y_n and x_n as exposed in the "Evaluation of the scores" section.
##
##########################
# FUNCTIONS
#
library(cglasso) # last cglasso'library version: download.file("https://github.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/raw/main/02%20-%20Simulations/SimulationStudy2/CGLASSO/cglasso_2.0.5.tar.gz", destfile = "cglasso_2.0.5.tar.gz")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/datajcggm.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/jcglasso.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/jcggm.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/gof.R")
source("https://raw.githubusercontent.com/gianluca-sottile/Hematopoiesis-network-inference-from-RT-qPCR-data/main/01%20-%20RCode/to_graph.R")
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


###############
# Import the data
load("IAGOS smooth data.RData") 

y.names <- c("O3", "NO", "H2O", "CO"); x.names <- "temp"
p <- length( y.names); q <- length(x.names)
N <- length(DF)

selected_T <- seq(1,max(df$alt), by = 50 )
N_T <- length(selected_T )
new_data <- data.frame(x = selected_T)

## # SIGMA_j MATRIX : N_t X N_T -dimensional covariance matrices for each variable:  
sigma_j <- vector( mode = "list", length = p+q)
for(j in seq_len(p+q) ){
  sigma_ji_num <- sigma_ji_den <- matrix(0, N_T,N_T)
  for(i in seq_len(N)){
    A <- weight[j, i, ] * t(YX.star[j, i, ] - YX.mean[j, ])
    sigma_ji_num <- sigma_ji_num + t(A) %*% A 
    sigma_ji_den <- sigma_ji_den + weight[j, i, ]%*% t(weight[j, i, ])
  }
  sigma_j[[j]] <- sigma_ji_num / sigma_ji_den
  print(j)
}

## # H - matrice H trece-class operato:
H <- matrix(0, N_T,N_T)
for (j in seq_len(p+q)) H = H + (1/(p+q))*sigma_j[[j]]

###  EIGEN():
eigen_dec <- eigen(H, symmetric = T)
# 99% variability
L <- which(cumsum(eigen_dec$values[ which(eigen_dec$values > 0)])/sum(eigen_dec$values[ which(eigen_dec$values > 0)]) > 0.99)[1]
L
phi <- eigen_dec$vectors[, c(1:L)]

### SCORES (centered Yjn - Yj.mean):
gamma_chi <- array(NA, dim = c( p+q, N, L ), dimnames = list( c(y.names, x.names), paste0("N_", seq(1, N) ), paste0("L_", seq(1, L) ) ) )
for( j in seq_len( p+q ) ){
  for( i in seq_len( N) ){
    if(sum(weight[c(y.names, x.names)[j], i, ]) != 0 ){
      gamma_chi[c(y.names, x.names)[j],i, ] <- solve(t(phi) %*% diag( weight[c(y.names, x.names)[j], i, ] ) %*% phi)  %*% t(phi) %*% diag( weight[c(y.names, x.names)[j], i, ] ) %*%( YX.star[c(y.names, x.names)[j], i, ] - YX.mean[c(y.names, x.names)[j],])
    } else gamma_chi[c(y.names, x.names)[j],i, ] <- 0
  }}

### EXPLAINED VARIABILITY at each element of the expansion:
n.matrix.list <- vector( mode = "list", length = L)
sd.matrix <- matrix(NA, p+q, L)
colnames( sd.matrix ) <- paste("l_", seq(1, L, by = 1))
rownames( sd.matrix ) <- c(y.names, x.names)
for( l in seq_len(L)){
  n.matrix.list[[l]] <- matrix(NA, N, p+q) # Scores of all units for all variables, for the expansion element 'l'
  for (n in seq_len(N)) {
    n.matrix.list[[l]][n,] <- gamma_chi[,n,l] # Scores of all variables for unit 'n' and element 'l'
  }
  sd.matrix[,l] <- apply(n.matrix.list[[l]], 2, sd)
}

### fGGRMs estimator

weights <- rep(n, L) / (n * L)
rho.min.ratio <- 0
rho.max.ratio <- 1
nrho <- 101
perc.rho.seq <- seq(from = rho.max.ratio , to= rho.min.ratio , length.out = nrho)
lmbd.min.ratio <- 0
lmbd.max.ratio <- 1
nlambda <- 101
perc.lmbd.seq <- seq(from= lmbd.max.ratio , to= lmbd.min.ratio , length.out = nlambda)

x <- y <- data <- eps.Y <- B.stim <- out <- vector(mode = "list", length = L)
l <- 1
for( l in seq_len(L)) {
  x[[l]] <- matrix(NA, N,q)
  colnames(x[[l]]) <- x.names
  y[[l]] <- matrix(NA, N,p)
  colnames(y[[l]]) <- y.names
  for (n in 1:N) {
    x[[l]][n,] <- gamma_chi[ x.names ,n, l]
    y[[l]][n,] <- gamma_chi[ y.names ,n, l]
  }
  y[[l]] <- t(apply(y[[l]], 1, function(n) (n / apply(y[[l]], 2, sd)) ))
  x[[l]] <- apply(x[[l]], 1, function(n) (n / apply(x[[l]], 2, sd)) )
  
  data <- datacggm( Y = y[[l]], X = x[[l]] )
  out[[l]] <- cglasso(.~. , data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  eps.Y[[l]] <- matrix(NA, n, q)
  B.stim[[l]] <- matrix(NA, p, q+1)# array(NA, dim = c(q, p, nrho) )
  eps.Y[[l]] <- y[[l]]- out[[l]]$mu[,,2,1]
  B.stim[[l]] <- t( out[[l]]$B[-1, , 2, ])
}
data.list <- datajcggm(Y = y, X = x)

out1 <- jcglasso(data = datajcggm(Y = eps.Y), nrho = 1, alpha1 = 0)
rho_max <- out1$rho

rho1 <- perc.rho.seq * rho_max
out.tht <- jcglasso(data = datajcggm(Y = eps.Y), rho = rho1, alpha1 = 0, nu = 0, alpha3 = 0)
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

niter <- 3
rho_estim <- lmbd_estim <- matrix(NA, 3, niter, dimnames = list( c("KLCV", "AIC", "eBIC"), niter = paste0("sim_", seq(1, niter)) ))
rho_estim[ "KLCV", 1] <- out.tht$rho[ id.min.klcv ]
rho_estim[ "AIC", 1] <- out.tht$rho[ id.min.aic ]
rho_estim[ "eBIC", 1] <- out.tht$rho[ id.min.ebic ]


iter <- 1
for ( iter in seq_len(niter-1)) {
  ### lambda KLCV
  lambda.klcv <- lambda.iteration( data = data.list , rho = rho_estim["KLCV" , iter], L = L, p = p,
                                   perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  
  b.stim <- vector( mode = "list", length = L)
  for(l in seq_len(L)){
    b.stim[[l]] <- array(NA, dim = c(p,1, length( perc.lmbd.seq) ))
    for (j in seq_along(perc.lmbd.seq)) {
      b.stim[[l]][,1,j] <- lambda.klcv$b.stim[[1]][,j] 
    }}
  
  id.min.klcv <- which.min( KLCV.gauss2(tht.stimYX = lambda.klcv$tht.stimYX, b.stim = b.stim, data = lambda.klcv$data, L = L, whole.tht = F ))
  lmbd_estim[ "KLCV", iter ] <- lambda.klcv$lambda[ id.min.klcv ]
  
  rm(lambda.klcv , id.min.klcv )
  ### lambda AIC
  lambda.aic <- lambda.iteration( data = data.list , rho = rho_estim[ "AIC", iter], L = L, p = p,
                                  perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  id.min.aic <- which.min( AIC.jcglasso(lambda.aic$out.jcglasso)$value_gof )
  lmbd_estim[ "AIC", iter ] <- lambda.aic$lambda[ id.min.aic ] 
  
  rm(lambda.aic , id.min.aic )
  
  ### lambda eBIC
  lambda.ebic <- lambda.iteration( data = data.list , rho = rho_estim["eBIC", iter], L = L,  p = p, 
                                   perc.lmbd.seq = perc.lmbd.seq, weights = weights)
  id.min.ebic <- which.min(BIC.jcglasso(lambda.ebic$out.jcglasso, g = .5, type = "FD")$value_gof )
  lmbd_estim["eBIC", iter ] <- lambda.ebic$lambda[id.min.ebic]
  
  rm(lambda.ebic , id.min.ebic )
  
  #
  #
  
  #### rho KLCV
  rho.klcv <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim[ "KLCV", iter], L = L, perc.rho.seq = perc.rho.seq)
  
  b.stim <- vector( mode = "list", length = L)
  for(l in seq_len(L)){
    b.stim[[l]] <- array(NA, dim = c(p,1, length( perc.rho.seq) ))
    for (j in seq_along(perc.rho.seq)) {
      b.stim[[l]][,1,j] <- rho.klcv$b.stim[[1]][,j] 
    }}
  
  id.min.klcv <- which.min( KLCV.gauss2( tht.stimYX = rho.klcv$tht.stim, b.stim = b.stim, data = rho.klcv$data, L = L, whole.tht = F ) )
  rho_estim["KLCV", iter + 1] <- rho.klcv$rho[id.min.klcv ]
  
  rm(rho.klcv , id.min.klcv)
  
  #### rho AIC
  rho.aic <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim["AIC", iter], L = L, perc.rho.seq = perc.rho.seq)
  id.min.aic <- which.min(AIC(rho.aic$out)$value_gof)
  rho_estim["AIC", iter + 1] <- rho.aic$rho[id.min.aic ]
  
  rm(rho.aic , id.min.aic)
  
  
  #### rho eBIC
  rho.ebic <- rho.iteration( y = y, x = x, n = n , p = p, q = q, lambda = lmbd_estim["eBIC", iter], L = L, perc.rho.seq = perc.rho.seq)
  id.min.ebic <- which.min(BIC(rho.ebic$out, g = .5, type = "FD")$value_gof)
  rho_estim[ "eBIC", iter + 1] <- rho.ebic$rho[id.min.ebic]
  
  rm(rho.ebic , id.min.ebic )
  
  iter <- iter + 1
}

### lambda KLCV
lambda.klcv <- lambda.iteration( data = data.list , rho = rho_estim["KLCV", niter], L = L, p = p,
                                 perc.lmbd.seq = perc.lmbd.seq, weights = weights)

b.stim <- vector( mode = "list", length = L)
for(l in seq_len(L)){
  b.stim[[l]] <- array(NA, dim = c(p,1, length( perc.lmbd.seq) ))
  for (j in seq_along(perc.lmbd.seq)) {
    b.stim[[l]][,1,j] <- lambda.klcv$b.stim[[1]][,j] 
  }}

id.min.klcv <- which.min( KLCV.gauss2(tht.stimYX = lambda.klcv$tht.stimYX, b.stim = b.stim, data = lambda.klcv$data, L = L, whole.tht = F ))
lmbd_estim["KLCV", niter] <- lambda.klcv$lambda[id.min.klcv]

KLCV.tht <- matrix(0, p,p)
KLCV.B <- matrix(0, q,p)
b.stim.plot <- list()
for (l in seq_len(L)) {
  KLCV.tht <- KLCV.tht + lambda.klcv$tht.stimY[[l]][,,id.min.klcv]
  KLCV.B <- KLCV.B + b.stim[[l]][,,id.min.klcv]
  b.stim.plot[[l]] <- b.stim[[l]][,,id.min.klcv]
}

rm(lambda.klcv , id.min.klcv, b.stim)

### lambda AIC
lambda.aic <- lambda.iteration( data = data.list , rho = rho_estim["AIC", niter], L = L, p = p,
                                perc.lmbd.seq = perc.lmbd.seq, weights = weights)

b.stim <- vector( mode = "list", length = L)
for(l in seq_len(L)){
  b.stim[[l]] <- array(NA, dim = c(p,1, length( perc.lmbd.seq) ))
  for (j in seq_along(perc.lmbd.seq)) {
    b.stim[[l]][,1,j] <- lambda.aic$b.stim[[1]][,j] 
  }}

id.min.aic <- which.min( AIC.jcglasso(lambda.aic$out)$value_gof )
lmbd_estim["AIC", niter] <- lambda.aic$lambda[id.min.aic]

AIC.tht <- matrix(0, p,p)
AIC.B <- matrix(0, q,p)
for (l in seq_len(L)) {
  AIC.tht <- AIC.tht + lambda.aic$tht.stimY[[l]][,,id.min.aic]
  AIC.B <- AIC.B + b.stim[[l]][,,id.min.aic]
  b.stim[[l]][,,id.min.aic]
}

rm(lambda.aic , id.min.aic, b.stim)

### lambda eBIC
lambda.ebic <- lambda.iteration( data = data.list , rho = rho_estim["eBIC", niter], L = L, p = p,
                                 perc.lmbd.seq = perc.lmbd.seq, weights = weights)

b.stim <- vector( mode = "list", length = L)
for(l in seq_len(L)){
  b.stim[[l]] <- array(NA, dim = c(p,1, length( perc.lmbd.seq) ))
  for (j in seq_along(perc.lmbd.seq)) {
    b.stim[[l]][,1,j] <- lambda.ebic$b.stim[[1]][,j] 
  }}

id.min.ebic <- which.min(BIC.jcglasso(lambda.ebic$out, g = .5, type = "FD")$value_gof )
lmbd_estim["eBIC", niter] <- lambda.ebic$lambda[id.min.ebic]

eBIC.tht <- matrix(0, p,p)
eBIC.B <- matrix(0, q,p)
for (l in seq_len(L)) {
  eBIC.tht <- eBIC.tht + lambda.ebic$tht.stimY[[l]][,,id.min.ebic]
  eBIC.B <- eBIC.B + b.stim[[l]][,,id.min.ebic]
}

rm(lambda.ebic , b.stim)

#############
# To visualize the results
library(igraph)

colnames(KLCV.tht) <- rownames(KLCV.tht) <- y.names
KLCV.tht[which(KLCV.tht != 0) ] <- 1
net.tht.KLCV <- graph_from_adjacency_matrix(KLCV.tht , mode = "lower", diag = F)
V(net.tht.KLCV )$color <- "white"
V(net.tht.KLCV )$size <- 34
V(net.tht.KLCV )$label.color <- "black"
V(net.tht.KLCV )$label.cex <- 1.5
E(net.tht.KLCV )$color <- "black"
V(net.tht.KLCV)$label <-  expression(O[3],  NO,  H2O, CO )
plot(net.tht.KLCV,layout=layout.circle ,asp = 1, main = "")
KLCV.B

colnames(AIC.tht) <- rownames(AIC.tht) <- y.names
AIC.tht[which(AIC.tht != 0) ] <- 1
net.tht.AIC <- graph_from_adjacency_matrix(AIC.tht , mode = "lower", diag = F)
V(net.tht.AIC )$color <- "white"
V(net.tht.AIC )$size <- 34
V(net.tht.AIC )$label.color <- "black"
V(net.tht.AIC )$label.cex <- 1.4
E(net.tht.AIC )$color <- "black"
V(net.tht.AIC)$label <-  expression(O[3],  NO,  H2O, CO)
plot(net.tht.AIC,layout=layout.circle ,asp = 1, main = "AIC")


colnames(eBIC.tht) <- rownames(eBIC.tht) <- y.names
eBIC.tht[which(eBIC.tht != 0) ] <- 1
net.tht.eBIC <- graph_from_adjacency_matrix(eBIC.tht , mode = "lower", diag = F)
V(net.tht.eBIC )$color <- "white"
V(net.tht.eBIC )$size <- 34
V(net.tht.eBIC )$label.color <- "black"
V(net.tht.eBIC )$label.cex <- 1.4
E(net.tht.eBIC )$color <- "black"
V(net.tht.eBIC)$label <-  expression(O[3],  NO,  H2O, CO)
plot(net.tht.eBIC,layout=layout.circle ,asp = 1, main = "eBIC")

colnames(KLCV.B) <- colnames(AIC.B) <- colnames(eBIC.B) <- y.names
as.numeric(KLCV.B != 0)
as.numeric(AIC.B != 0)
as.numeric(eBIC.B != 0)
KLCV.B

B.hat <- matrix(0, N_T, N_T)
B.hat.l <- 0
for( l in seq_len(L)){
  B.hat.l <- b.stim.plot[[l]][3] * (phi[,l] %*% t(phi[,l]))
  B.hat <- B.hat + B.hat.l
}

max(B.hat); min(B.hat)
labels <- new_data$x

## PLOT estimated B(t,t'):
par(mfrow = c(1,1), cex.lab = 5)
par( cex.lab = 1)
heatmap( B.hat, labRow = labels, labCol = labels, Rowv = NA, Colv = NA, ylab = "temp - altitude", xlab = "altitude", margins = c(1.5,1.5))
