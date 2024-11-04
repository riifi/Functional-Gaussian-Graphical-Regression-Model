load("IAGOS DF.Rdata")

df <- bind_rows(DF)

y.names <- c("O3", "NO", "H2O", "CO"); x.names <- "temp"
p <- length( y.names); q <- length(x.names)

DF1 <- DF
N <- length(DF1)

DF <- vector( mode = "list", length = N)
id.rm <- id.no.var <- c()
for( n in seq_len(N) ){
  if ( dim(DF1[[n]])[2] < 9) {
    DF[[n]] <- NULL
    id.no.var <- append(id.no.var, n)
  }
  else{
    a <- DF1[[n]]
    a[DF1[[n]] == -999999.9] <- NULL
    DF[[n]] <- na.omit(a)
    DF[[n]] <-  DF[[n]][ which( DF[[n]]$alt < 13000 
                                & DF[[n]]$CO < mean(DF[[n]]$CO)+ 3* sd(DF[[n]]$CO)
                                & DF[[n]]$H2O < mean(DF[[n]]$H2O)+ 3* sd(DF[[n]]$H2O)
                                & DF[[n]]$NO < mean(DF[[n]]$NO)+ 3* sd(DF[[n]]$NO)
                                & DF[[n]]$O3 < mean(DF[[n]]$O3)+ 3 * sd(DF[[n]]$O3))
                         ,]
    if ( dim(DF[[n]])[1] <1 ) id.rm <- append( id.rm, n)
    print( paste0(n, " - ", dim(DF[[n]])[1] ))
  }
}
DF2 <- DF[!sapply(DF, is.null)]
N <- length(DF2)
df <- bind_rows(DF2)
sum(is.na(df)); sum( df == -999999.9)

a <- c(0)
for( n in 1:N) print( paste0(n, " - ", dim(DF2[[n]]) ))
for( n in 1:N) if( dim(DF2[[n]])[1] == 0 ) a <- append(a, n)
a
DF <- DF2[-a]
N <- length(DF)

df <- bind_rows(DF)

rm(nc_files, DF1, DF2, file_names)


###################
N <- length(DF)

selected_T <- seq(1,max(df$alt), by = 50 )
N_T <- length(selected_T )
new_data <- data.frame(x = selected_T)

weight <- Y.star <- array(NA, dim = c( "var"= p + q, "unit" = N, "Tpoint" = N_T ), dimnames = list( c(y.names, x.names), paste0("N_", seq(1, N) ), paste0("T_", seq(1, N_T) ) ) )
check.weight <- weight <- YX.star <- array(NA, dim = c( "var"= p + q , "unit" = N, "Tpoint" = N_T ), dimnames = list( c(y.names, x.names), paste0("N_", seq(1, N) ), paste0("T_", seq(1, N_T) ) ) )
YX.mean <- matrix(0, p +q, N_T); rownames(YX.mean) <- c(y.names,x.names); colnames(YX.mean) <- paste0("T_", seq(1, N_T) ) 
for( j in seq_len(p+q) ){
  ## calculate the observations at delected points stored in new_data:
  for( i in seq_len(N)){
    y <- DF[[i]][, c(y.names, x.names)[j] ]
    x <- DF[[i]][,"alt" ]
    data <- data.frame( y = y, x = x)
    model_i <- gam(y ~ s(x), data = data)
    predictions <- predict(model_i, newdata = new_data, type = "response", se.fit = TRUE)
    w_ji <- 1/((predictions$se.fit + 1)^2)
    weight[c(y.names, x.names)[j] , i, ] <- w_ji
    YX.star[c(y.names, x.names)[j] , i, ] <- predictions$fit
    YX.mean[c(y.names, x.names)[j] , ] <- YX.mean[c(y.names, x.names)[j] , ] + weight[c(y.names, x.names)[j] , i, ] * YX.star[c(y.names, x.names)[j], i, ]
  }
  YX.mean[c(y.names, x.names)[j] , ] <- YX.mean[c(y.names, x.names)[j] , ] / apply( weight[c(y.names, x.names)[j] , , ], 2, sum)
}

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
heatmap( B.hat, labRow = labels, labCol = labels, Rowv = NA, Colv = NA, ylab = "temp - altitude", xlab = "altitude", margins = c(1.5,1.5))#,  label.cex = 4)
