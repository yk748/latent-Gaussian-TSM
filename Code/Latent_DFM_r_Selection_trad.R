#################################################################################
# Function for computing information criteria and edge distribution for selecting r
#################################################################################
Latent_DFM_r_Selection_trad <- function(d,TT,p,Cov_Z,identy_opt){
  
  # ------------------------------------------------- #
  # set up for three ICs and ED
  r_grid <- seq(1,10)
  IC <- array(NA,dim=c(length(r_grid),3))
  iter_max <- 4
  r_hat_delta <- vector("numeric",length=iter_max)
  
  # ------------------------------------------------- #
  # main loop for ICs
  for (r in r_grid){
    
    Estim_r <- Latent_DFM_Estim(r,p,Cov_Z,identy_opt,shift=TRUE)
    Cov_eta_r <- diag(Estim_r$Cov_eps,d)
    V_term <- log(norm(Cov_eta_r,type="F")^2/(d*TT))
    C_dT <- min(sqrt(d),sqrt(TT))
    
    IC[r,1] <- V_term + r*(d+TT)/(d*TT)*log(d*TT/(d+TT))
    IC[r,2] <- V_term + r*(d+TT)/(d*TT)*log(C_dT^2)
    IC[r,3] <- V_term + r*log(C_dT^2)/C_dT^2
  }
  min_idx_IC <- apply(IC,2,which.min)
  r_IC1 <- min_idx_IC[1]
  r_IC2 <- min_idx_IC[2]
  r_IC3 <- min_idx_IC[3]
  
  # ------------------------------------------------- #
  # set up for ED
  iter_max <- 10
  r_hat_delta <- rep(NA,iter_max)
  
  Cov_Z0 <- Cov_Z[,,(p+1)]
  if (sum(eigen(Cov_Z0)$values < 0) > 0){
    eps <- 10e-8
    shift_Cov_Z0 <- Cov_Z0 + (1+eps)*abs(eigen(Cov_Z0)$values[d])*diag(1,d)
  }else{
    shift_Cov_Z0 <- Cov_Z0
  }
  scaler <- diag(shift_Cov_Z0)
  Cov_Z0 <- shift_Cov_Z0 / (sqrt(scaler) %*% t(sqrt(scaler)))
  eig_val <- eigen(Cov_Z0)$values
  
  # ------------------------------------------------- #
  # main loop for ED
  r_max <- 10
  j <- r_max +1
  r_delta <- vector("numeric",length=iter_max)
  for (iter in 1:iter_max){
    
    y <- eig_val[j:(j+4)]
    x <- mapply(x = c(0:4),FUN = function(x) (j-1+x)^(2/3))
    beta_hat <- solve(t(x) %*% x) %*% t(x) %*% y
    delta <- 2*abs(beta_hat)
    
    cond <- which(-diff(eig_val) >= as.numeric(delta))
    if ( identical(cond,integer(0)) ){
      break
    }
    
    r_hat_delta[iter] <- max(which(cond <= r_max))
    if ( identical(r_hat_delta[iter], integer(0)) ){
      r_hat_delta[iter] <- 0
    }
    
    j <- r_delta[iter] + 1
  }
  r_ED <- r_delta[iter]
  
  
  output <- list()
  output$r_IC1 <- r_IC1
  output$r_IC2 <- r_IC2
  output$r_IC3 <- r_IC3
  output$r_ED <- r_ED
  return(output)
}