#################################################################################
# Function for computing information criteria for selecting p
#################################################################################
Latent_DFM_p_Selection_trad <- function(d,TT,r,X_t,Link,identy_opt){
  
  # ------------------------------------------------- #
  # set up for four ICs
  p_max <- 5
  p_grid <- seq(1,p_max)
  IC <- array(NA,dim=c(length(p_grid),4))
  
  # ------------------------------------------------- #
  # Compute Cov_X and implied Cov_Z with p_max
  Cov_X_p_max <- Latent_Gauss_Cov(d,TT,p_max,X_t)
  Cov_Z_p_max <- Latent_Gauss_invLink(d,p,Cov_X_p_max,Link)
  
  # ------------------------------------------------- #
  # main loop for ICs
  for (p in p_grid){
    
    # ------------------------------------------------- #
    # Estimate model parameters with Cov_Z within p
    Estim_p <- Latent_DFM_Estim(r,p,Cov_Z_p_max[,,c((p_max+1-p):(p_max+1+p))],
                                identy_opt,shift=TRUE)
    Cov_eta_p <- diag(Estim_p$Cov_eps,d)
    V_term <- log(norm(Cov_eta_p,type="F")^2/(d*TT))
    
    # AIC
    IC[p,1] <- V_term + 2/TT*p*r^2
    
    # SC
    IC[p,2] <- V_term + 2*log(log(TT))/TT*p*r^2
    
    # HQ
    IC[p,3] <- V_term + log(TT)/TT*p*r^2
    
    # FPE
    IC[p,4] <- V_term + 2*r*(r*p+1)/TT
  }
  min_idx_IC <- apply(IC,2,which.min)
  p_IC1 <- min_idx_IC[1]
  p_IC2 <- min_idx_IC[2]
  p_IC3 <- min_idx_IC[3]
  p_IC4 <- min_idx_IC[4]
  
  output <- list()
  output$p_IC1 <- p_IC1
  output$p_IC2 <- p_IC2
  output$p_IC3 <- p_IC3
  output$p_IC4 <- p_IC3
  return(output)
}
