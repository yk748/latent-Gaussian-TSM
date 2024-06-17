#################################################################################
# Main function for generating discrete-valued observations and latent dynamic factor model
#################################################################################
Latent_DFM_Model <- function(p,r,d,TT,dist,dist_opt=NULL,model,identy_opt){
  
  # ------------------------------------------------- #
  # Set up
  Psi_star <- model$Psi
  Cov_eta_star <- model$Cov_eta
  Cov_eps_star <- model$Cov_eps
  Lambda_star <- model$Lambda
  Param <- model$Param

  # ------------------------------------------------- #
  # Compute Sigma_Y(0) and compute scalers
  if (p == 1){
    Cov_U_star <- Cov_eta_star
    C_star <- matrix(Psi_star,r,r)
    
  }else{
    Cov_U_star <- matrix(0,p*r,p*r)
    Cov_U_star[1:r,1:r] <- Cov_eta_star
    C_star <- as.matrix(companion_form_phi(Psi_star,r,p))
    
  }
  
  diff <- 2
  order <- 0
  Cov_G_star <- Cov_U_star
  while( diff > 1){
    order <- order + 1
    if ( r == 1 ){
      update_mat <- C_star^order %*% Cov_U_star %*% t(C_star^order)
    }else{
      update_mat <- C_star^order %*% Cov_U_star %*% t(C_star^order) 
    }
    
    Cov_G_star <- Cov_G_star + update_mat
    if ( r == 1 ){
      diff <- update_mat
    }else{
      diff <- norm(update_mat,type="o")
    }
  }
  
  Cov_Y0_star <- Cov_G_star[1:r,1:r]
  common_term <- Lambda_star %*% Cov_Y0_star %*% t(Lambda_star)
  Cov_Z0_star <- common_term + Cov_eps_star
  scaler <- sqrt(diag(Cov_Z0_star)) 

  # ------------------------------------------------- #
  # Scale other parameters
  Cov_eps <- Cov_eps_star / (scaler %*% t(scaler))
  
  if (identy_opt == 1){
    eigen <- eigen(common_term/ (scaler %*% t(scaler)))
    Ur <- eigen$vectors[,1:r]
    scale_f <- sqrt(eigen$values[1:r] / diag(Cov_Y0_star))
    Lambda <- Ur %*% diag(scale_f,r)
    
    YK <- Cov_G_star - C_star %*% Cov_G_star %*% t(C_star)
    Cov_eta <- YK[1:r,1:r]
    Psi <- Psi_star
    
  }else if (identy_opt == 2){
    Lambda <- (Lambda_star / scaler) %*% solve((Lambda_star / scaler)[1:r,1:r])
    Cov_eta <- ((Lambda_star / scaler)[1:r,1:r]) %*% Cov_eta_star %*% t(((Lambda_star / scaler)[1:r,1:r]))
    Psi <- array(NA,dim=dim(Psi_star))
      if (p == 1){
        Psi[,,1] <- ((Lambda_star / scaler)[1:r,1:r]) %*% Psi_star[,,1] %*% solve(((Lambda_star / scaler)[1:r,1:r]))
      }else{
        for (h in 1:p){
          Psi[,,h] <- ((Lambda_star / scaler)[1:r,1:r]) %*% Psi_star[,,h] %*% solve(((Lambda_star / scaler)[1:r,1:r]))
        }
      }
  }  
  
  # ------------------------------------------------- #
  # generate sample paths of Y_t and Z_t
  Burn <- TT
  Y_t <- array(NA, dim=c(r,TT+Burn))
  Z_t <- array(NA, dim=c(d,TT+Burn))
  
  eps_t <- t(mvtnorm::rmvnorm(TT+Burn, rep(0,d), Cov_eps, method="eigen"))
  eta_t <- t(mvtnorm::rmvnorm(TT+Burn, rep(0,r), Cov_eta, method="eigen"))
  
  for (t in 1:(TT+Burn) ){
    if (t <= p){
      
      Y_t[,t] <- eta_t[,t]
      
    }else{
      if ( r == 1){
        Y_t[,t] <- sum(sapply(c(1:p),function(x){sum(Psi[,,x] %*% Y_t[,(t-x)])})) +eta_t[,t]
      }else{
        Y_t[,t] <- rowSums(sapply(c(1:p),function(x){rowSums(Psi[,,x] %*% Y_t[,(t-x)])})) +eta_t[,t] 
      }
    }
    
    Z_t[,t] <- Lambda %*% as.matrix(Y_t[,t]) + eps_t[,t] 
  }
  Y_t <- Y_t[,-c(1:Burn)]
  Z_t <- Z_t[,-c(1:Burn)]
  
  # ------------------------------------------------- #
  # compute thresholds
  if ( dist == "Bern"){
    Prob <- unlist(Param)
    
    threshold <- lapply(Prob, function(x){qnorm(1-x)})
    min_val <- 0
    
  } else if( dist == "Pois") { 
    Theta <- unlist(Param)
    
    range_max <- sapply(Theta, function(x){which(ppois(0:100,x,log=FALSE)>=1-10^(-8))[1]})
    threshold <- mapply(x=1:d,function(x){qnorm(ppois(0:range_max[x],Theta[x],log=FALSE))})
    min_val <- 0
    
  } else if( dist == "multinom"){
    
    num_category <- dist_opt
    threshold <- lapply(X=1:d,FUN=function(X)qnorm(cumsum(Param[[X]][1:(num_category-1)]),0,1))
    min_val <- 1
    
  } else if (dist == "negbin"){
    Prob <- mapply(x=1:d,function(x)Param[[x]][1])
    success <- mapply(x=1:d,function(x)Param[[x]][2])
    
    range_max <- sapply(Prob, function(x){which(pnbinom(0:500,size=success,x,log.p = FALSE)>=1-10^(-8))[1]})
    threshold <- mapply(x=1:d,function(x){qnorm(pnbinom(0:range_max[x],size=success[x],Prob[x],log=FALSE))})
    min_val <- 0
    
  }
  
  # ------------------------------------------------- #
  # generate sample paths of X_t
  X_t <- array(NA, dim=c(d,TT))
  if (dist == "Bern"){
    for (i in 1:d){
      thrs <- threshold[[i]]
      X_t[i,] <- ( Z_t[i,] > thrs )*1
      
    } 
    
  } else{
    if (dist == "multinom"){
      min_val <- 1
    }else{ # dist == "Pois" | dist == "negbin"
      min_val <- 0
    }
    
    for (i in 1:d){
      aug_thrs <- c(-Inf, threshold[[i]], Inf)
      for (t in 1:TT){
        X_t[i,t] <- max(cumsum((Z_t[i,t] > aug_thrs))) -1 + min_val
      }
    }
  }
  

  output <- list()
  output$Lambda <- Lambda
  output$Cov_eps <- Cov_eps
  output$Psi <- Psi
  output$Cov_eta <- Cov_eta
  output$threshold <- threshold
  output$Y_t <- Y_t
  output$Z_t <- Z_t
  output$X_t <- X_t
  return(output)
}


#################################################################################
# Function that converts VAR transition matrices into an augmented VAR(1) transition matrix
#################################################################################
companion_form_phi = function(Psi,d,p){
  
  comp_Psi <- array(0,dim=c((p*d),(p*d)))
  if (p == 1){
    comp_Psi = Psi[,,1]
  }
  else{
    for (i in 1:p){
      for (j in 1:p){
        if(i == 1){
          comp_Psi[1:d,((j-1)*d+1):(j*d)] <- Psi[,,j]
        }
        else if ( (i-1) == j & i != 1 ){
          comp_Psi[((i-1)*d+1):(i*d),((j-1)*d+1):(j*d)] <- diag(1,d)
        }
      }
    }
  }
  
  return(comp_Psi)
}
