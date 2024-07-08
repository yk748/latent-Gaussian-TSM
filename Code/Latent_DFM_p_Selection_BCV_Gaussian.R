#################################################################################
# Function for block cross-validation methods for selecting p with Z_t
#################################################################################
Latent_DFM_p_Selection_BCV_Gaussian <- function(d,TT,r,Z_t,identy_opt,fold=4){
  
  # ------------------------------------------------- #
  # set up for PC
  p_max <- 5
  p_grid <- seq(1,p_max)
  MSE_PC <- array(NA,dim=c(length(p_grid),fold))
  
  # ------------------------------------------------- #
  # Data splitting index
  idx <- c(1,sapply(X=1:(fold-1),FUN=function(X) 1 + X*floor(TT/fold) ),(TT+1))
  
  # ------------------------------------------------- #
  # main loop together
  for (m in 1:fold){
    # ------------------------------------------------- #
    # Data splitting
    test_idx <- idx[m]:(idx[m+1]-1)
    train_idx <- c(1:TT)[-c(test_idx)]
    Z_t_train <- Z_t[,train_idx]
    Z_t_test <- Z_t[,test_idx]
    
    TT_train <- dim(Z_t_train)[2]
    TT_test <- dim(Z_t_test)[2]
    
    # ------------------------------------------------- #
    # Construct Cov_Z_hat from train dataset
    Cov_Z_train_p_max <- Latent_Gauss_Cov(d,TT_train,p_max,Z_t_train)
    
    # ------------------------------------------------- #
    # Construct implied Cov_Z_hat from test dataset
    Cov_Z_test_p_max <- Latent_Gauss_Cov(d,TT_test,p_max,Z_t_test)
    
    for (p in p_grid){
      
      # ------------------------------------------------- #
      # Estimation with train dataset
      Estim_train_p <- Latent_DFM_Estim(r,p,Cov_Z_train_p_max[,,c((p_max+1-p):(p_max+1+p))],identy_opt,shift=FALSE)
      
      Psi_train_p <- array(NA,dim=c(p*r,r))
      for (h in 1:p){
        Psi_train_p[((h-1)*r+1):(h*r),] <- t(Estim_train_p$Psi[,,h])
      }
      vec_Psi_train_p <- as.vector(Psi_train_p)
      
      # Estimation with test dataset
      Estim_test_p <- Latent_DFM_Estim(r,p,Cov_Z_test_p_max[,,c((p_max+1-p):(p_max+1+p))],identy_opt,shift=FALSE)
      
      Gamma_test_Y0 <- array(NA,dim=c(p*r,r))
      for (h in 1:p){
        Gamma_test_Y0[((h-1)*r+1):(h*r),] <- t(Estim_test_p$Cov_Y[,,(p+1+h)])
      }
      gamma_test_Y <- as.vector(Gamma_test_Y0)
      
      Gamma_test_Yp <- array(NA,dim=c(p*r,p*r))
      for (i in 1:p){
        for (j in 1:p){
          if ( i < j ){
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Estim_test_p$Cov_Y[,,(p+1-abs(i-j))])
          }else if( i == j ){
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- Estim_test_p$Cov_Y[,,(p+1)]
          }else{ # i > j
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Estim_test_p$Cov_Y[,,(p+1+abs(i-j))])
          }
        }
      }
      Gamma_test_Y <- bdiag(replicate(r,Gamma_test_Yp,simplify=FALSE))
      
      # ------------------------------------------------- #
      # ------------------------------------------------- #
      # Compute MSE
      MSE_PC[p,m] <-  as.numeric(-2*t(vec_Psi_train_p)%*%gamma_test_Y + t(vec_Psi_train_p)%*%Gamma_test_Y%*%vec_Psi_train_p)
    }
  }
  
  MSE_p <- rowSums(MSE_PC)/fold
  idx_min <- which.min(MSE_p)
  p_hat <- p_grid[idx_min]
  
  output <- list()
  output$p_PC <- p_hat
  output$MSE_PC <- MSE_p
  return(output)
}