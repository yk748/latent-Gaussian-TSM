#################################################################################
# Function for the latent dynamic factor model-based forecasting
#################################################################################
Latent_DFM_Forecast <- function(model,p,r,d,TT,X_t,H,dist,dist_opt){
  
  # ------------------------------------------------- #
  # Set up from models
  Lambda <- model$Lambda
  Cov_eps = model$Cov_eps
  Psi <- companion_form_phi(model$Psi,r,p)
  Cov_eta <- model$Cov_eta
  
  if (is.null(model$threshold)){
    threshold <- NULL
  }else{
    threshold <- model$threshold
  }
  Param <- model$Param
  
  # ------------------------------------------------- #
  # Set up for Kalman recursions
  window_sample <- 5
  N <- 100
  window <- min(TT,window_sample)
  Y_hat <- array(NA,dim=c(p*r,window,N))
  Q_hat <- array(NA,dim=c(p*r,p*r,window))
  Z_hat <- array(NA,dim=c(d,window,N))
  R_hat <- array(NA,dim=c(d,d,window))
  
  Y_tilde <- array(NA,dim=c(p*r,(window+1),N))
  Q_tilde <- array(NA,dim=c(p*r,p*r,(window+1)))
  
  Psi_n <- diag(1,p*r)
  Q_init <- Cov_eta
  for (n in 1:10){
    Psi_n <- Psi %*% Psi_n
    Q_init <- Q_init + Psi_n %*% Cov_eta %*% t(Psi_n)
  }
  
  Q_tilde[,,1] <- Q_init
  Y_tilde[,1,] <- rmvnorm(N,mean=rep(0,p*r),sigma=Q_init,method="eigen")
  
  # ------------------------------------------------- #
  # Set up for particle filtering
  N <- 100 # number of particles
  Z_t <- array(NA,dim=c(d,window ,N))
  w_t <- array(NA,dim=c((window+1),N))
  w_t[1,] <- 1
  w_t[1,] <- w_t[1,]/sum(w_t[1,])
  alpha_t <- array(NA,dim=c(window ,N))
  ESS <- array(NA,dim=window)
  
  # ------------------------------------------------- #
  # compute thresholds & upper and lower bounds
  if (is.null(threshold)){
    if ( dist == "Bern"){
      Prob <- unlist(Param)
      threshold <- lapply(Prob, function(x){qnorm(1-x)})
      min_val <- 0
      
    } else if( dist == "multinom") { 
      num_category <- dist_opt
      threshold <- lapply(X=1:d,FUN=function(X)qnorm(cumsum(Param[[X]][1:(num_category-1)]),0,1))
      min_val <- 1
      
    } else if( dist == "Pois"){
      Theta <- unlist(Param)
      range_max <- sapply(Theta, function(x){which(ppois(0:100,x,log=FALSE)>=1-10^(-8))[1]})
      threshold <- mapply(x=1:d,function(x){qnorm(ppois(0:range_max[x],Theta[x],log=FALSE))})
      min_val <- 0
      
    } else if (dist == "negbin"){
      Prob <- mapply(x=1:d,function(x)Param_hat[[x]][1])
      success <- mapply(x=1:d,function(x)Param_hat[[x]][2])
      range_max <- sapply(Prob, function(x){which(pnbinom(0:500,size=success,x,log.p = FALSE)>=1-10^(-8))[1]})
      threshold <- mapply(x=1:d,function(x){qnorm(pnbinom(0:range_max[x],size=success[x],Prob[x],log=FALSE))})
      min_val <- 0
      
    } 
  }else{
    if ( dist == "Bern"){
      min_val <- 0
    }else if(dist == "multinom"){
      min_val <- 1
    }else if( dist == "Pois"){
      min_val <- 0
    }else if(dist == "negbin"){
      min_val <- 0
    }
  }
  
  lower_bdd <- array(NA,dim=c(d,window))
  upper_bdd <- array(NA,dim=c(d,window))
  M <- Inf
  for (i in 1:d){
    aug_thrs <- c(-M,threshold[[i]],M)
    val <- X_t[i,(TT-window+1):(TT)] - min_val
    lower_bdd[i,] <- aug_thrs[val+1]
    upper_bdd[i,] <- aug_thrs[val+2]
  }
  
  # ------------------------------------------------- #
  # main loop
  for (t in 1:window){
    # Forecast step
    Q_hat[,,t] <- Psi %*% Q_tilde[,,t] %*% t(Psi) + Cov_eta
    R_hat[,,t] <- Lambda %*% Q_hat[,,t] %*% t(Lambda) + Cov_eps
    
    for (k in 1:N){
      Y_hat[,t,k] <- Psi %*% Y_tilde[,t,k]
      Z_hat[,t,k] <- Lambda %*% Y_hat[,t,k]
    }
    
    # ------------------------------------------------- #
    # SIS step
    for (k in 1:N){
      Z_t[,t,k] <- TruncatedNormal::rtmvnorm(n=1,mu=Z_hat[,t,k],sigma=R_hat[,,t],
                                             lb=lower_bdd[,t],ub=upper_bdd[,t])
      
      alpha_t[t,k] <- TruncatedNormal::pmvnorm(mu=Z_hat[,t,k],sigma=R_hat[,,t],
                                               lb=lower_bdd[,t],ub=upper_bdd[,t])[1]
      
      w_t[(t+1),k] <- w_t[t,k]*alpha_t[t,k]
    }
    w_t[(t+1),] <- w_t[(t+1),]/sum(w_t[(t+1),])
    
    # ------------------------------------------------- #
    # resampling
    ESS[t] <- 1/sum((w_t[(t+1),])^2)
    if (ESS[t] < N/2){
      W <- c(0, cumsum(w_t[(t+1),]))
      U <- (0:(N-1))/N + runif(1,min=0,max=1/N)
      idx_sys <- sapply(seq_along(U), function(i) (sum(W[i] <= U & U < W[i + 1])))
      new_idx <- rep(1:N,idx_sys)
      
      Y_hat[,t,] <- Y_hat[,t,new_idx]
      Z_hat[,t,] <- Z_hat[,t,new_idx]
      Z_t[,t,] <- Z_t[,t,new_idx]
      w_t[(t+1),] <- 1/N
    }
    
    # ------------------------------------------------- #
    # Update step
    R_hat_inv <- solve(R_hat[,,t])
    K_gain <- Q_hat[,,t]%*%t(Lambda)%*%R_hat_inv
    Q_tilde[,,(t+1)] <- (diag(1,r) - K_gain%*%Lambda)%*%Q_hat[,,t]
    
    for (k in 1:N){
      Y_tilde[,(t+1),k] <- Y_hat[,t,k] + K_gain%*%(Z_t[,t,k] - Z_hat[,t,k])
    }
  }
  Y_tilde <- Y_tilde[,2:(window+1),]
  Q_tilde <- Q_tilde[,,2:(window+1)]
  w_t <- w_t[2:(window+1),]
  
  # ------------------------------------------------- #
  # Prediction after the observation period
  Y_hat_h <- array(NA,dim=c(r,H,N))
  Q_hat_h <- array(NA,dim=c(r,r,H))
  Z_hat_h <- array(NA,dim=c(d,H,N))
  R_hat_h <- array(NA,dim=c(d,d,H))
  
  Y_tilde_h <- Y_tilde[,window,]
  Q_tilde_h <- Q_tilde[,,window]
  w_T <- w_t[window,]
  
  for (h in 1:H){
    for (k in 1:N){
      if (h == 1){
        Y_hat_h[,h,k] <- Psi %*% Y_tilde_h[,k]
        Q_hat_h[,,h] <- Psi %*% Q_tilde_h  %*% t(Psi) + Cov_eta
        
      }else{
        Y_hat_h[,h,k] <- Psi %*% Y_hat_h[,(h-1),k]
        Q_hat_h[,,h] <- Psi %*% Q_hat_h[,,(h-1)]  %*% t(Psi) + Cov_eta
      }
      Z_hat_h[,h,k] <- Lambda %*% Y_hat_h[,h,k]
      R_hat_h[,,h] <- Lambda %*% Q_hat_h[,,h] %*% t(Lambda) + Cov_eps
    }
  }
  
  # ------------------------------------------------- #
  # Create discrete-valued observations
  X_t_h <- array(NA,dim=c(d,H))
  
  for (h in 1:H){
    for (i in 1:d){
      D_V_i <- array(NA,dim=c(length(threshold[[i]])+1,N))
      
      for (k in 1:N){
        for (l in 1:(length(threshold[[i]])+1)){
          if (l == 1){
            lower_bdd <- -Inf
            upper_bdd <- threshold[[i]][l]
            
          }else if (l > 1 & l <= length(threshold[[i]])){
            lower_bdd <- threshold[[i]][l-1]
            upper_bdd <- threshold[[i]][l]
            
          }else{
            lower_bdd <- threshold[[i]][l-1]
            upper_bdd <- Inf
            
          }
          D_V_i[l,k] <- c( pnorm(upper_bdd,Z_hat_h[i,h,k],R_hat_h[i,i,h]) 
                            - pnorm(lower_bdd,Z_hat_h[i,h,k],R_hat_h[i,i,h]) )
        }
      }
      X_t_h[i,h] <- which.max(w_T %*% t(D_V_i)) - 1 + min_val
    }
  }
  
  output <- list()
  output$Y_hat <- Y_hat
  output$Y_tilde <- Y_tilde
  output$Q_hat <- Q_hat
  output$Q_tilde <- Q_tilde
  output$Z_hat <- Z_hat
  output$R_hat <- R_hat
  output$threshold <- threshold
  output$Z_t <- Z_t
  output$w_t <- w_t
  output$alpha_t <- alpha_t
  output$ESS <- ESS
  output$X_t_h <- X_t_h
  output$Z_hat_h <- Z_hat_h 
  output$Y_hat_h <- Y_hat_h
  return(output)
}
