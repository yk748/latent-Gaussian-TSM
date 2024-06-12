#################################################################################
# Function for computing autocovariance matrices of discrete-valued observations
#################################################################################
Latent_Gauss_Cov <- function(d,TT,p,X_t){
  
  # ------------------------------------------------- #
  # Estimate Cov_X
  Cov_X <- array(NA,dim=c(d,d,(2*p+1)))
  for (i in 1:d){
    for (j in 1:d){
      for (h in 1:(2*p+1)){
        if (h < p){
          idx_front <- (h+1):TT
          idx_back <- 1:(TT-h)
          
          xt_i <- X_t[i,idx_front]
          if (sum(xt_i==0)==length(idx_front)){
            xt_i[1] <- 1e-8
          }
          xt_j <- X_t[j,idx_back]
          if (sum(xt_j==0)==length(idx_back)){
            xt_j[1] <- 1e-8
          }
          
          Cov_X[i,j,h] <- cor(xt_i,xt_j)*(TT-h)/TT
          
        }else if ( h == p+1 ){
          Cov_X[i,j,h] <- cor(X_t[i,],X_t[j,])
          
        }else{ # h > p
          idx_front <- 1:(TT-h)
          idx_back <- (h+1):TT
          
          xt_i <- X_t[i,idx_front]
          if (sum(xt_i==0)==length(idx_front)){
            xt_i[1] <- 1e-8
          }
          xt_j <- X_t[j,idx_back]
          if (sum(xt_j==0)==length(idx_back)){
            xt_j[1] <- 1e-8
          }
          
          Cov_X[i,j,h] <- cor(xt_i,xt_j)*(TT-h+p+1)/TT
          
        }
      }
    }
  }
  # 
  # # ------------------------------------------------- #
  # # Verify if any entries in Cov_X is NA
  # for (h in 1:(2*p+1)){
  #  if ( h == p+1){
  #    next
  #    
  #  }else{
  #    for (j in 1:d){
  #      sum_front <- sum(is.na(Cov_X[,j,h]))
  #      sum_back <- sum(is.na(Cov_X[,j,(2*p+2-h)]))
  #      if(sum_front !=0 | sum_back !=0){
  #        if (sum_front!=0){
  #          Cov_X[,j,h] <- Cov_X[,j,(2*p+2-h)]
  #        }else{ # sum_back !=0
  #          Cov_X[,j,(2*p+2-h)] <- Cov_X[,j,h]
  #        }
  #      }
  #      
  #    }
  #  }
  #   
  # }
  
  return(Cov_X)
}
