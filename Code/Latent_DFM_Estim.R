#################################################################################
# Function for estimating the latent dynamic factor model
#################################################################################
Latent_DFM_Estim <- function(r,p,Cov_Z,identy_opt,shift){
  
  # ------------------------------------------------- #
  # PCA setup
  Cov_Z0 <- Cov_Z[,,(p+1)]
  if (shift==TRUE){
    if (sum(eigen(Cov_Z0)$values < 0) > 0){
      eps <- 10e-8
      shift_Cov_Z0 <- Cov_Z0 + (1+eps)*abs(eigen(Cov_Z0)$values[d])*diag(1,d)
    }else{
      shift_Cov_Z0 <- Cov_Z0
    }
    scaler <- diag(shift_Cov_Z0)
    Cov_Z0 <- shift_Cov_Z0 / (sqrt(scaler) %*% t(sqrt(scaler)))
  }
  
  eigen_decomp <- eigen(Cov_Z[,,(p+1)])
  Ur <- eigen_decomp$vector[,1:r]
  Dr <- diag(eigen_decomp$values[1:r],r)
  
  # ------------------------------------------------- #
  # Estimate the loadings matrix
  if (identy_opt == 1){
    Xi <- Ur %*% sqrt(Dr)
    Lambda <- Xi
    
  }else if (identy_opt == 2){
    Xi <- Ur %*% sqrt(Dr)
    if (r == 1){
      Lambda <- Xi / Xi[1,1]
    }else{
      Lambda <- Xi %*% solve(Xi[1:r,1:r])
    }
    
  }
  
  # ------------------------------------------------- #
  # Estimate the covariance matrix of the noise
  Cov_eps <- diag(Cov_Z0 - Xi %*% t(Xi))

  # ------------------------------------------------- #
  # Estimate the covariance matrices of the factor models
  Cov_Y <- array(NA,dim=c(r,r,(2*p+1)))
  Q_inv <- solve(t(Lambda)%*%Lambda)
  for (h in 1:(2*p+1)){
    if (h == (p+1)){
      if (identy_opt == 1){
        Cov_Y[,,(p+1)] <- diag(1,r)
      }else if(identy_opt == 2){
        Cov_Y[,,(p+1)] <- Xi[1:r,1:r] %*% t(Xi[1:r,1:r])
      }
      
    }else{
      if (identy_opt == 1){
        Cov_Y[,,h] <- solve(sqrt(Dr)) %*% t(Ur) %*% Cov_Z[,,h] %*% Ur %*% solve(sqrt(Dr))
      }else if(identy_opt == 2){
        Cov_Y[,,h] <- Q_inv %*% ( t(Lambda) %*% Cov_Z[,,h] %*% Lambda ) %*% Q_inv  
      }
      
    }
  }  
 
  # ------------------------------------------------- #
  # Estimate the VAR transition matrices of the factor series
  Gamma_Y0 <- array(NA,dim=c(p*r,r))
  for (h in 1:p){
    Gamma_Y0[((h-1)*r+1):(h*r),] <- t(Cov_Y[,,(p+1+h)])
  }
  Gamma_Yp <- array(NA,dim=c(p*r,p*r))
  for (i in 1:p){
    for (j in 1:p){
      if ( i < j ){
        Gamma_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y[,,(p+1-abs(i-j))])
      }else if( i == j ){
        Gamma_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- Cov_Y[,,(p+1)]
      }else{ # i > j
        Gamma_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y[,,(p+1+abs(i-j))])
      }
    }
  }
  YW_sol <- solve(Gamma_Yp + diag(1e-8,p*r)) %*% Gamma_Y0
  
  Psi <- array(NA,dim=c(r,r,p))
  for (h in 1:p){
    Psi[,,h] <- t(YW_sol[((h-1)*r+1):(h*r),])
  }
  
  # ------------------------------------------------- #
  # Estimate the covariance matrices of the factor series
  tmp_hat <- array(0,dim=c(r,r))
  for (h in 1:p){
    tmp_hat <- tmp_hat + Psi[,,h] %*% t(Cov_Y[,,(p+1+h)])
  }
  Cov_eta <- Cov_Y[,,(p+1)] - tmp_hat
  
  
  output <- list()
  output$Lambda <- Lambda
  output$Cov_eps <- Cov_eps
  output$Cov_Y <- Cov_Y 
  output$Psi <- Psi
  output$Cov_eta <- Cov_eta
  return(output)
}
