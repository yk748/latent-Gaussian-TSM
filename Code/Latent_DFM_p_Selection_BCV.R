#################################################################################
# Function for block cross-validation methods for selecting p
#################################################################################
Latent_DFM_p_Selection_BCV <- function(d,TT,r,X_t,Link,identy_opt,fold=4){
  
  # ------------------------------------------------- #
  # set up for PC
  p_grid <- seq(1,5); p_max <- 5
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
    X_t_train <- X_t[,train_idx]
    X_t_test <- X_t[,test_idx]
    
    TT_train <- dim(X_t_train)[2]
    TT_test <- dim(X_t_test)[2]
    
    # ------------------------------------------------- #
    # Construct implied Cov_Z_hat from train dataset
    Cov_X_train_p_max <- Latent_Gauss_Cov(d,TT_train,p_max,X_t_train)
    nan <- sum(is.nan(Cov_X_train_p_max))!=0
    na <- sum(is.na(Cov_X_train_p_max))!=0
    if (nan | na){
      Cov_X_p_max <- Latent_Gauss_Cov(d,TT,p_max,X_t)
      
      for (h in 1:(2*p+1)){
        if (nan){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.nan(Cov_X_train_p_max[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.nan(Cov_X_train_p_max[,j,h])))!=0)  
        }else if(na){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.na(Cov_X_train_p_max[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.na(Cov_X_train_p_max[,j,h])))!=0)
        }
        
        for (i in mis_idx_i){
          Cov_X_train_p_max[i,,h] <- Cov_X_p_max[i,,h]
        }
        for (j in mis_idx_j){
          Cov_X_train_p_max[,j,h] <- Cov_X_p_max[,j,h]
        }
      }
      
    }
    Cov_Z_train_p_max <- Latent_Gauss_invLink(d,p_max,Cov_X_train_p_max,Link)
    
    # ------------------------------------------------- #
    # Construct implied Cov_Z_hat from test dataset
    Cov_X_test_p_max <- Latent_Gauss_Cov(d,TT_test,p_max,X_t_test)
    nan <- sum(is.nan(Cov_X_test_p_max))!=0
    na <- sum(is.na(Cov_X_test_p_max))!=0
    if (nan | na){
      Cov_X_p_max <- Latent_Gauss_Cov(d,TT,p_max,X_t)
      
      for (h in 1:(2*p+1)){
        if (nan){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.nan(Cov_X_test_p_max[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.nan(Cov_X_test_p_max[,j,h])))!=0)  
        }else if(na){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.na(Cov_X_test_p_max[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.na(Cov_X_test_p_max[,j,h])))!=0)
        }
        
        for (i in mis_idx_i){
          Cov_X_test_p_max[i,,h] <- Cov_X_p_max[i,,h]
        }
        for (j in mis_idx_j){
          Cov_X_test_p_max[,j,h] <- Cov_X_p_max[,j,h]
        }
      }
    }
    Cov_Z_test_p_max <- Latent_Gauss_invLink(d,p_max,Cov_X_test_p_max,Link)
    
    
    for (p in p_grid){
      
      # ------------------------------------------------- #
      # Estimation with train dataset
      Estim_train_p <- Latent_DFM_Estim(r,p,Cov_Z_train_p_max[,,c((p_max+1-p):(p_max+1+p))],identy_opt,shift=TRUE)
      
      Psi_train_p <- array(NA,dim=c(p*r,r))
      for (h in 1:p){
        Psi_train_p[((h-1)*r+1):(h*r),] <- t(Estim_train_p$Psi[,,h])
      }
      vec_Psi_train_p <- as.vector(Psi_train_p)
      
      # Estimation with test dataset
      Estim_test_p <- Latent_DFM_Estim(r,p,Cov_Z_test_p_max[,,c((p_max+1-p):(p_max+1+p))],identy_opt,shift=TRUE)
      
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