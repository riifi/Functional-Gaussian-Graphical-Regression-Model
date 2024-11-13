##########
## TO RUN IF YOU WANT TO GET THE REALIZATIONS y_n and x_n AS INDICATED IN SECTION 
## "EVALUATION OF THE SCORES" FROM THE FILES AIRCRAFT FILES DOWNLOADED FROM https://iagos.aeris-data.fr/download/
##########

library(ncdf4)
library(dplyr)
library(mgcv)

# Download L2 files for the year 2020 from https://iagos.aeris-data.fr/download/
# Each file is a flight, and each flight is a statistical unit

folder_path <- "Desktop/IAGOS/2020 L2" # Set YOUR folder path
file_names <- list.files(path = folder_path, pattern = "\\.nc4$", full.names = TRUE)

# Open all the files .nc4
nc_files <- lapply(file_names, nc_open)

# Select the files with the variables we are interested in analysing and store them in a list
DF <- vector(mode = "list", length = length(nc_files))
for( j in seq_along(DF)){
  try.alt <- try(ncvar_get( nc_files[[j]], "gps_alt_AC") )
  try.NO <- try(ncvar_get( nc_files[[j]], "NO_P2b") )
  if( !inherits(try.alt, "try-error") & !inherits(try.NO, "try-error")){
    DF[[j]] <- as.data.frame(cbind( "lat" = ncvar_get( nc_files[[j]], "lat"),
                                    "lon" = ncvar_get( nc_files[[j]], "lon"),
                                    "alt" = ncvar_get( nc_files[[j]], "gps_alt_AC"),
                                    "temp" = ncvar_get( nc_files[[j]], "air_temp_AC"),
                                    "CO" = ncvar_get( nc_files[[j]], "CO_P1"),
                                    "O3" = ncvar_get( nc_files[[j]], "O3_P1"),
                                    "NO" = ncvar_get( nc_files[[j]], "NO_P2b") ,
                                    "H2O"= ncvar_get( nc_files[[j]], "H2O_gas_P1")
    ))
    DF[[j]]$time <- rep( nc_files[[j]]$dim$UTC_time$units, length(ncvar_get( nc_files[[j]], "baro_alt_AC")))
  }
  else print(j)
}
DF <- DF[!sapply(DF, is.null)]
df <- bind_rows(DF)

y.names <- c("O3", "NO", "H2O", "CO"); x.names <- "temp"
p <- length( y.names); q <- length(x.names)

DF1 <- DF
N <- length(DF1)

# Remove flights with missing variables, remove outliers and missing data indicated by the value -999999.9
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

# Remove flights without observation after outlier and NA cleaning
a <- c(0)
for( n in 1:N) print( paste0(n, " - ", dim(DF2[[n]]) ))
for( n in 1:N) if( dim(DF2[[n]])[1] == 0 ) a <- append(a, n)
a
DF <- DF2[-a]

# Total number of remaining flights
N <- length(DF)

df <- bind_rows(DF)

rm(nc_files, DF1, DF2, file_names, folder_path,id.no.var, id.rm)

N <- length(DF)

# Fix the grid
selected_T <- seq(1,max(df$alt), by = 50 )
N_T <- length(selected_T )
new_data <- data.frame(x = selected_T)


# Generalized additive model function to have the smoothed data for the selected grid
weight <- Y.star <- array(NA, dim = c( "var"= p + q, "unit" = N, "Tpoint" = N_T ), dimnames = list( c(y.names, x.names), 
                                                                                                    paste0("N_", seq(1, N) ), paste0("T_", seq(1, N_T) ) ) )
check.weight <- weight <- YX.star <- array(NA, dim = c( "var"= p + q , "unit" = N, "Tpoint" = N_T ), dimnames = list( c(y.names, x.names), paste0("N_", seq(1, N) ), paste0("T_", seq(1, N_T) ) ) )
YX.mean <- matrix(0, p +q, N_T); rownames(YX.mean) <- c(y.names,x.names); colnames(YX.mean) <- paste0("T_", seq(1, N_T) ) 
for( j in seq_len(p+q) ){
  ## Observations of variables at selected points
  for( i in seq_len(N)){
    y <- DF[[i]][, c(y.names, x.names)[j] ]
    x <- DF[[i]][,"alt" ] 
    data <- data.frame( y = y, x = x)
    model_i <- gam(y ~ s(x), data = data) 
    predictions <- predict(model_i, newdata = new_data, type = "response", se.fit = TRUE) # plot(new_data[,1], predictions$fit, cex = 0.3)
    w_ji <- 1/((predictions$se.fit + 1)^2) 
    weight[c(y.names, x.names)[j] , i, ] <- w_ji
    YX.star[c(y.names, x.names)[j] , i, ] <- predictions$fit
    
    YX.mean[c(y.names, x.names)[j] , ] <- YX.mean[c(y.names, x.names)[j] , ] + weight[c(y.names, x.names)[j] , i, ] * YX.star[c(y.names, x.names)[j], i, ]
    
  }
  YX.mean[c(y.names, x.names)[j] , ] <- YX.mean[c(y.names, x.names)[j] , ] / apply( weight[c(y.names, x.names)[j] , , ], 2, sum)
}

#save( df, DF, YX.mean, YX.star, weight,  file = "IAGOS smooth data.RData")

# Continue with "fGGRM IAGOS DATA.R" file
