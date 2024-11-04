library(cglasso)
library(mgcv)
library(ggplot2)
library(maps)

## To download the data:
# https://www.kaggle.com/datasets/adityaramachandran27/world-air-quality-index-by-city-and-coordinates/data
aqi <- read.csv("AQI and Lat Long of Countries.csv",header=TRUE, sep=",")

head(aqi)
str(aqi)

min.points.coord <- 100 # we select the units which present at least 100 measurement points
y.index <- c(3,5,7,9)
y.names <- colnames(aqi)[y.index]
p <- length(y.index)
x.index <- c(11)
x.names <- colnames(aqi)[x.index]
q <- length(x.index)
n.var <- length(y.index)+ length(x.index)

lng_window <- 3
int <- seq( min(aqi$lng),max(aqi$lng),lng_window )

summary(aqi)

aqi_ordered <- aqi[ order( aqi[, "lat"]), ]
aqi_ordered <- aqi_ordered[, c(y.index, x.index, which( colnames(aqi_ordered) == "lat" | colnames(aqi_ordered) == "lng"))]
aqi1 <- aqi
aqi <- aqi_ordered

## Define the statistical units according to longitude
yx <- NULL
cnt <- 1
ind <- NULL
points.n <- NULL
coord <- NULL
lat_N <- c()
for (i in 1:(length(int)-1)){
  print( sum(( aqi$lng > int[i] ) & ( aqi$lng < int[i+1] )))
  if (sum((aqi$lng > int[i] ) & (aqi$lng < int[i+1])) > min.points.coord){
    print(i)
    ind[cnt] <- i
    row.ind <- which( aqi$lng > int[i] & aqi$lng < int[i+1])
    yx[[cnt]] <- aqi[ row.ind , ]
    points.i <- sum( length(row.ind) )
    coord <- rbind(coord, aqi[ row.ind,  c("lng", "lat")])
    points.n[cnt] <- points.i
    cnt <- cnt+1
  }
}
length(ind)
N <- length(yx)

yx.unlist <- yx[[1]]
for( n in seq_len(N-1)) yx.unlist <- rbind(yx.unlist,yx[[n+1]])
#### PLOTS - MAPS
# world_map <- map_data("world")
# 
# dim(yx.unlist)
# ggplot() +
#   geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#bfbebd", color = "#bfbebd") +
#   geom_point(data = yx.unlist, aes(x = lng, y = lat), color = "black", size = 0.2) +
#   labs(x = "Longitude", y = "Latitude") +
#   theme_minimal()+
#   theme( axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
#          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14) )
#
# ggplot() +
#   geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#bfbebd", color = "#bfbebd") +
#   geom_point(data = aqi, aes(x = lng, y = lat), color = "#160863", size = 0.3, alpha = 0.5) +
#   geom_point(data = yx.new.unlist, aes(x = lng, y = lat), color = "black", size = 0.3) +
#   labs(x = "Longitude", y = "Latitude") +
#   theme_minimal()

## seleziono i punti del dominio su cui calcolare le fitted functions
# plot(density(aqi$lat))
# hist(aqi$lat, nclass = 100)
# min.dens.lat <- which(aqi$lat > - 40)[1]
# max.dens.lat <- which(aqi$lat >  70)[1] -1
# selected_T <- aqi$lat[seq(min.dens.lat, max.dens.lat, by = 10 )]

## Grid for the predictions
selected_T <- sort(yx.unlist$lat)[seq(1, dim(yx.unlist)[1], by = 10 )][-1]
N_T <- length(selected_T )

##  PLOT - selected grid
# plot.data1 <- expand.grid( selected_T, int[ind] + lng_window/2 )
# colnames(plot.data1)
# ggplot() +
#   geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#bfbebd", color = "#bfbebd") +
#   geom_point(data = plot.data1, aes(x = Var2, y = Var1), color = "black", size = 0.000001) +
#   labs(x = "Longitude", y = "Latitude") +
#   theme_minimal()+
#   theme( axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
#          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14) )

weight <- Y.star <- array(NA, dim = c( "var"=p + q, "unit" = N, "Tpoint" = N_T ), dimnames = list( c(y.names, x.names), 
                                                                                                   paste0("N_", seq(1, N) ), paste0("T_", seq(1, N_T) ) ) )
Y.mean <- matrix(NA, p+q, N_T)
new_data <- data.frame(x = selected_T)
for( j in seq_len(p+q) ){
  y <- aqi[, c(y.names, x.names)[j] ]
  x <- aqi[, "lat" ] # plot(x,y, cex = 0.3)
  model <- gam( y ~ s(x))
  ## Once an estimate of the variable's value is obtained considering [-70,70], T = N = total number of observations, this represents the average curve at the selected points in the domain:
  Y.mean[j, ] <- predict(model, newdata = new_data, type = "response")  # plot(new_data[,1],Y.mean[j, ], cex = .3) 
  ## Calculate the observations of the variables at the selected points:
  for( i in seq_len(N) ){
    y <- yx[[i]][, c(y.names, x.names)[j] ]
    x <- yx[[i]][,"lat" ] # plot(x,y, cex = 0.3)
    data <- data.frame( y = y, x = x)
    model_i <- gam(y ~ s(x), data = data) # plot(model_i) 
    predictions <- predict(model_i, newdata = new_data, type = "response", se.fit = TRUE) # plot(new_data[,1], predictions$fit, cex = 0.3)
    ## Predictions represent the observed values without error (smoothing is applied) at the selected points. 
    ## Calculate the weights for each selected point in T for the unit:
    w_ji <- 1/(predictions$se.fit^2)
    if( any(predictions$se.fit == 0) ){
      w_ji[which( predictions$se.fit == 0 )] <- 0
    }
    weight[j, i, ] <- w_ji
    Y.star[j, i, ] <- predictions$fit
  }
}

## # SIGMA_j MATRIX : N_t X N_T -dimensional covariance matrices for each variable:  
sigma_j <- vector( mode = "list", length = p+q)
for(j in seq_len(p+q) ){
  sigma_ji_num <- sigma_ji_den <- matrix(0, N_T,N_T)
  for(i in seq_len(N)){
    A <- weight[j, i, ] * t(Y.star[j, i, ] - Y.mean[j, ])
    sigma_ji_num <- sigma_ji_num + t(A) %*% A 
    sigma_ji_den <- sigma_ji_den + weight[j, i, ]%*% t(weight[j, i, ])
  }
  sigma_j[[j]] <- sigma_ji_num / sigma_ji_den
  print(j)
}

## # H - matrice H trece-class operato:
H <- matrix(0, N_T,N_T)
for (j in seq_len(p+q)) H = H + (1/(p+q))*sigma_j[[j]]

######  EIGEN():
eigen_dec <- eigen(H, symmetric = T)
# 99% variability
L <- which(cumsum(eigen_dec$values[ which(eigen_dec$values > 0)])/sum(eigen_dec$values[ which(eigen_dec$values > 0)]) > 0.99)[1]
phi <- eigen_dec$vectors[, c(1:L)]

## # SCORES (centered Yjn - Yj.mean):
gamma_chi <- array(NA, dim = c( p+q, N, L ), dimnames = list( c(y.names, x.names), paste0("N_", seq(1, N) ), paste0("L_", seq(1, L) ) ) )
for( j in seq_len( p+q ) ){
  for( i in seq_len( N) ){
    if(sum(weight[j, i, ]) != 0 ){
      gamma_chi[j,i, ] <- solve(t(phi) %*% diag( weight[j, i, ] ) %*% phi)  %*% t(phi) %*% diag( weight[j, i, ] ) %*%( Y.star[j, i, ] - Y.mean[j,])
    } else gamma_chi[j,i, ] <- 0
  }}


## # EXPLAINED VARIABILITY at each element of the expansion:
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
## PLOT s.d.:
# par(mfrow = c(1,1))
# par(mar = c(5, 5, 5, 5), cex.main = 2)
# matplot(t(sd.matrix), type = "l", lwd = 2.3, xlab = "Terms of the expansion", ylab = "s.d.", cex.axis = 1.3, cex.lab = 1.7)
# legend("topright", legend = c( expression(AQI, CO, O[3], NO[2]), "PM 2.5" ), lwd = 2.3,  col = 1:5, lty = 1:5, cex = 1.3)

weights <- rep(n, L) / (n * L)
rho.min.ratio <- 0
rho.max.ratio <- 1
nrho <- 31
perc.rho.seq <- seq(from = rho.max.ratio , to= rho.min.ratio , length.out = nrho)
lmbd.min.ratio <- 0
lmbd.max.ratio <- 1
nlambda <- 31
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

# plot(density( data$X[,1] ))
# plot(density( data$Y[,2]  ))
# plot(density(   data$Y[,3]  ))
# plot(density(   data$Y[,4]  ))

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

# Update the labels with chemical elements

colnames(KLCV.tht) <- rownames(KLCV.tht) <- c( "AQI", "CO", "O3", "NO2")
KLCV.tht[which(KLCV.tht != 0) ] <- 1
net.tht.KLCV <- graph_from_adjacency_matrix(KLCV.tht , mode = "lower", diag = F)
V(net.tht.KLCV )$color <- "white"
V(net.tht.KLCV )$size <- 34
V(net.tht.KLCV )$label.color <- "black"
V(net.tht.KLCV )$label.cex <- 2.5
E(net.tht.KLCV )$color <- "black"
V(net.tht.KLCV)$label <-  expression(AQI, CO, O[3], NO[2])
plot(net.tht.KLCV,layout=layout.circle ,asp = 1, main = "")
KLCV.B

colnames(AIC.tht) <- rownames(AIC.tht) <- c( "AQI", "CO", "O3", "NO2")
AIC.tht[which(AIC.tht != 0) ] <- 1
net.tht.AIC <- graph_from_adjacency_matrix(AIC.tht , mode = "lower", diag = F)
V(net.tht.AIC )$color <- "white"
V(net.tht.AIC )$size <- 34
V(net.tht.AIC )$label.color <- "black"
V(net.tht.AIC )$label.cex <- 1.4
E(net.tht.AIC )$color <- "black"
V(net.tht.AIC)$label <-  expression(AQI, CO, O[3], N[O[2]])
plot(net.tht.AIC,layout=layout.circle ,asp = 1, main = "AIC")


colnames(eBIC.tht) <- rownames(eBIC.tht) <- c( "AQI", "CO", "O3", "NO2")
eBIC.tht[which(eBIC.tht != 0) ] <- 1
net.tht.eBIC <- graph_from_adjacency_matrix(eBIC.tht , mode = "lower", diag = F)
V(net.tht.eBIC )$color <- "white"
V(net.tht.eBIC )$size <- 34
V(net.tht.eBIC )$label.color <- "black"
V(net.tht.eBIC )$label.cex <- 1.4
E(net.tht.eBIC )$color <- "black"
V(net.tht.eBIC)$label <-  expression(AQI, CO, O[3], N[O[2]])
plot(net.tht.eBIC,layout=layout.circle ,asp = 1, main = "eBIC")


#### PLOT - NETS
# par(mfrow = c(1,3))
# par(mar = c(25, 0, 6, 2), cex.main = 2)
# plot(net.tht.KLCV,layout=layout.circle ,asp = 1, main = "jKLCV")
# plot(net.tht.AIC,layout=layout.circle ,asp = 1, main = "AIC")
# plot(net.tht.eBIC,layout=layout.circle ,asp = 1, main = "eBIC")

colnames(KLCV.B) <- colnames(AIC.B) <- colnames(eBIC.B) <- c( "AQI", "CO", "O3", "NO2")
as.numeric(KLCV.B != 0)
as.numeric(AIC.B != 0)
as.numeric(eBIC.B != 0)
KLCV.B

B.hat <- matrix(0, N_T, N_T)
B.hat.l <- 0
for( l in seq_len(L)){
  B.hat.l <- b.stim.plot[[l]][1] * (phi[,l] %*% t(phi[,l]))
  B.hat <- B.hat + B.hat.l
}

max(B.hat); min(B.hat)

## PLOT estimated B(t,t'):
# par(mfrow = c(1,1), cex.lab = 5)
# par( cex.lab = 5)
# heatmap( B.hat, Rowv = NA, Colv = NA, labRow = NA, labCol = NA, ylab = "PM2.5 latitude", xlab = "AQI latitude", margins = c(1.5,1.5),  label.cex = 4)


