#################################################################################
# Function for block cross-validation methods for selecting p with Y_t
#################################################################################
Latent_DFM_p_Selection_BCV_Gaussian_lowdim <- function(TT,r,Y_t,fold=4){
  
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
    Y_t_train <- matrix(Y_t[,train_idx],nrow=r)
    Y_t_test <- matrix(Y_t[,test_idx],nrow=r)
    
    TT_train <- dim(Y_t_train)[2]
    TT_test <- dim(Y_t_test)[2]
    
    # ------------------------------------------------- #
    # Construct Cov_Y_hat from train dataset
    Cov_Y_train_p_max <- Latent_Gauss_Cov(r,TT_train,p_max,Y_t_train)
    
    # ------------------------------------------------- #
    # Construct implied Cov_Y_hat from test dataset
    Cov_Y_test_p_max <- Latent_Gauss_Cov(r,TT_test,p_max,Y_t_test)
    
    for (p in p_grid){
      
      # ------------------------------------------------- #
      # Estimate the VAR transition matrices of the factor series with train dataset
      Cov_Y_train_p <- array(Cov_Y_train_p_max[,,c((p_max+1-p):(p_max+1+p))],dim=c(r,r,(2*p+1)))
      
      Gamma_train_Y0 <- array(NA,dim=c(p*r,r))
      for (h in 1:p){
        Gamma_train_Y0[((h-1)*r+1):(h*r),] <- t(Cov_Y_train_p[,,(p+1+h)])
      }
      Gamma_train_Yp <- array(NA,dim=c(p*r,p*r))
      for (i in 1:p){
        for (j in 1:p){
          if ( i < j ){
            Gamma_train_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y_train_p[,,(p+1-abs(i-j))])
          }else if( i == j ){
            Gamma_train_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- Cov_Y_train_p[,,(p+1)]
          }else{ # i > j
            Gamma_train_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y_train_p[,,(p+1+abs(i-j))])
          }
        }
      }
      YW_sol_p <- solve(Gamma_train_Yp + diag(1e-8,p*r)) %*% Gamma_train_Y0
      
      Psi_train_p <- array(NA,dim=c(r,r,p))
      for (h in 1:p){
        Psi_train_p[,,h] <- t(YW_sol_p[((h-1)*r+1):(h*r),])
      }
      vec_Psi_train_p <- as.vector(Psi_train_p)

      # ------------------------------------------------- #
      # Estimation with test dataset
      Cov_Y_test_p <- array(Cov_Y_test_p_max[,,c((p_max+1-p):(p_max+1+p))],dim=c(r,r,(2*p+1)))
      
      Gamma_test_Y0 <- array(NA,dim=c(p*r,r))
      for (h in 1:p){
        Gamma_test_Y0[((h-1)*r+1):(h*r),] <- t(Cov_Y_test_p[,,(p+1+h)])
      }
      gamma_test_Y <- as.vector(Gamma_test_Y0)
      
      Gamma_test_Yp <- array(NA,dim=c(p*r,p*r))
      for (i in 1:p){
        for (j in 1:p){
          if ( i < j ){
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y_test_p[,,(p+1-abs(i-j))])
          }else if( i == j ){
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- Cov_Y_test_p[,,(p+1)]
          }else{ # i > j
            Gamma_test_Yp[((i-1)*r+1):(i*r),((j-1)*r+1):(j*r)] <- t(Cov_Y_test_p[,,(p+1+abs(i-j))])
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