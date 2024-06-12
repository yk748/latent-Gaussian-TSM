#################################################################################
# Function for block cross-validation methods for selecting r
#################################################################################
Latent_DFM_r_Selection_BCV <- function(d,TT,p,X_t,Link,identy_opt,fold=4){
  
  # ------------------------------------------------- #
  # set up for both PC and Fac
  r_grid <- seq(1,10)
  MSE_PC <- array(NA,dim=c(length(r_grid),fold))
  MSE_Fac <- array(NA,dim=c(length(r_grid),fold))
  
  # ------------------------------------------------- #
  # Data splitting index for both PC and Fac
  idx <- c(1,sapply(X=1:(fold-1),FUN=function(X) 1 + X*floor(TT/fold) ),(TT+1))
  
  # ------------------------------------------------- #
  # main loop together
  for (m in 1:fold){
    # ------------------------------------------------- #
    # Data splitting for both PC and Fac
    test_idx <- idx[m]:(idx[m+1]-1)
    train_idx <- c(1:TT)[-c(test_idx)]
    X_t_train <- X_t[,train_idx]
    X_t_test <- X_t[,test_idx]
    
    TT_train <- dim(X_t_train)[2]
    TT_test <- dim(X_t_test)[2]
    
    # ------------------------------------------------- #
    # Construct implied Cov_Z from train dataset for both PC and Fac
    Cov_X_train_r <- Latent_Gauss_Cov(d,TT_train,p,X_t_train)
    nan <- sum(is.nan(Cov_X_train_r))!=0
    na <- sum(is.na(Cov_X_train_r))!=0
    if (nan | na){
      Cov_X <- Latent_Gauss_Cov(d,TT,p,X_t)
      
      for (h in 1:(2*p+1)){
        if (nan){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.nan(Cov_X_train_r[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.nan(Cov_X_train_r[,j,h])))!=0)  
        }else if(na){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.na(Cov_X_train_r[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.na(Cov_X_train_r[,j,h])))!=0)
        }
        
        for (i in mis_idx_i){
          Cov_X_train_r[i,,h] <- Cov_X[i,,h]
        }
        for (j in mis_idx_j){
          Cov_X_train_r[,j,h] <- Cov_X[,j,h]
        }
      }
      
    }
    Cov_Z_train_r <- Latent_Gauss_invLink(d,p,Cov_X_train_r,Link)
    
    # ------------------------------------------------- #
    # Construct implied Gamma_Z_hat from test dataset for both PC and Fac
    Cov_X_test_r <- Latent_Gauss_Cov(d,TT_test,p,X_t_test)
    nan <- sum(is.nan(Cov_X_test_r))!=0
    na <- sum(is.na(Cov_X_test_r))!=0
    if (nan | na){
      Cov_X <- Latent_Gauss_Cov(d,TT,p,X_t)
      
      for (h in 1:(2*p+1)){
        if (nan){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.nan(Cov_X_test_r[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.nan(Cov_X_test_r[,j,h])))!=0)  
        }else if(na){
          mis_idx_i <- which(mapply(i=1:d,function(i)sum(is.na(Cov_X_test_r[i,,h])))!=0)
          mis_idx_j <- which(mapply(j=1:d,function(j)sum(is.na(Cov_X_test_r[,j,h])))!=0)
        }
        
        for (i in mis_idx_i){
          Cov_X_test_r[i,,h] <- Cov_X[i,,h]
        }
        for (j in mis_idx_j){
          Cov_X_test_r[,j,h] <- Cov_X[,j,h]
        }
      }
    }
    Cov_Z_test_r <- Latent_Gauss_invLink(d,p,Cov_X_test_r,Link)
    
    for (r in r_grid){
      # ------------------------------------------------- #
      # ------------------------------------------------- #
      # Estimation with train dataset for Fac
      Cov_Z0_r <- Cov_Z_train_r[,,(p+1)]
      eps <- 10e-8
      shift_Cov_Z0_r <- Cov_Z0_r + (1+eps)*abs(eigen(Cov_Z0_r)$values[d])*diag(1,d)
      scaler <- diag(shift_Cov_Z0_r)
      Cov_Z0_r <- shift_Cov_Z0_r / sqrt(scaler) %*% t(sqrt(scaler))
      
      psych_fit <- suppressMessages(psych::fa(r = Cov_Z0_r, nfactors = r, rotate = "oblimin", warnings = FALSE))
      if (r == 1){
        mi_cor <- psych_fit$loadings %*% t(psych_fit$loadings) + diag(diag(psych_fit$residual),d)
      }else{
        mi_cor <- psych_fit$loadings %*% psych_fit$Phi %*% t(psych_fit$loadings) + diag(diag(psych_fit$residual),d)
      }
      
      # ------------------------------------------------- #
      # ------------------------------------------------- #
      # Estimation with train dataset for PC
      Estim_train_r <- Latent_DFM_Estim(r,p,Cov_Z_train_r,identy_opt,shift=TRUE)
      Lambda_train_r <- Estim_train_r$Lambda
      Cov_Y0_train_r <- Estim_train_r$Cov_Y
      Cov_signal_train_r <- matrix(Lambda_train_r,d,r) %*% Cov_Y0_train_r[,,(p+1)] %*% t(matrix(Lambda_train_r,d,r))
      Cov_eta_train_r <- diag(Estim_train_r$Cov_eps,d)
      
      # ------------------------------------------------- #
      # ------------------------------------------------- #
      # Compute MSE for PC
      MSE_PC[r,m] <- norm(Cov_Z_test_r[,,(p+1)] - (Cov_signal_train_r + Cov_eta_train_r) ,type="F")^2
      
      # ------------------------------------------------- #
      # ------------------------------------------------- #
      # Compute MSE for Fac
      MSE_Fac[r,m] <- norm(Cov_Z_test_r[,,(p+1)] - mi_cor,"F")^2
    }
  }
  
  # Choose the minimizers for both PC and Fac
  MSE_Fac <- rowSums(MSE_Fac)/fold
  MSE_PC <- rowSums(MSE_PC)/fold
  idx_min_Fac <- which.min(MSE_Fac)
  idx_min_PC <- which.min(MSE_PC)
  r_Fac <- r_grid[idx_min_Fac]
  r_PC <- r_grid[idx_min_PC]
  
  output <- list()
  output$r_Fac <- r_Fac
  output$r_PC <- r_PC
  output$MSE_Fac <- MSE_Fac
  output$MSE_PC <- MSE_PC
  return(output)
}