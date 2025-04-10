#################################################################################
# Main function for computing pairwise link functions
#################################################################################
Latent_Gauss_Link <- function(d,TT,dist,X_t,dist_opt){
  
  # ------------------------------------------------- #
  # Estimate parameters of marginal distributions
  Param_hat <- list()
  for (i in 1:d){
    
    if (dist == "Bern"){
      Param_hat[[i]] <- sum(X_t[i,])/TT
      
    } else if (dist == "Pois"){ 
      Param_hat[[i]] <- sum(X_t[i,])/TT
      
    } else if (dist == "multinom"){
      num_category <- dist_opt
      Param_hat[[i]] <- mapply(x=1:num_category,function(x)sum(X_t[i,]==x)/TT)
      
    } else if (dist == "negbin"){
      size <- 3
      prob <- size / (size + sum(X_t[i,])/TT)
      # prob <- 1 - mean(X_t[i,])/var(X_t[i,])
      # size <- mean(X_t[i,])^2/(var(X_t[i,]) - mean(X_t[i,]))
      Param_hat[[i]] <- c(prob,size)
      
    }
  }
  
  # ------------------------------------------------- #
  # Estimate coefficients of link function
  K <- 100
  ell_ij_hat <- array(NA,dim=c(K,d,d));
  for (i in 1:d){
    param_list <- list()
    
    if (dist == "Bern"){
      param_list[[1]] <- Param_hat[[i]]
      
    }else if (dist == "Pois"){
      param_list[[1]] <- Param_hat[[i]]
      
    }else if (dist == "multinom"){
      param_list[[1]] <- Param_hat[[i]]
      param_list[[2]] <- num_category
      
    }else if (dist == "negbin"){
      param_list[[1]] <- Param_hat[[i]][1]
      param_list[[2]] <- Param_hat[[i]][2]
      
    }
    
    for (j in 1:d){
      if (j >= i){
        if (dist == "Bern"){
          param_list[[2]] <- Param_hat[[j]]
          
        }else if (dist == "Pois"){
          param_list[[2]] <- Param_hat[[j]]
          
        }else if (dist == "multinom"){
          param_list[[3]] <- Param_hat[[j]]
          param_list[[4]] <- num_category
          
        }else if (dist == "negbin"){
          param_list[[3]] <- Param_hat[[j]][1]
          param_list[[4]] <- Param_hat[[j]][2]
          
        }
        ell_ij_hat[,i,j] <- l_ij_func(param_list,dist,K)
        
      }else{
        ell_ij_hat[,i,j] <- ell_ij_hat[,j,i]
        
      }
      
    }
  }
  
  output <- list()
  output$Param <- Param_hat
  output$Link <- ell_ij_hat
  return(output)
}

#################################################################################
# Function for computing Hermite coefficients of the link functions
#################################################################################
l_ij_func <- function(param_list,dist,K){
  
  # ------------------------------------------------- #
  # ith entry:
  N <- 10000
  if (dist == "Bern"){
    P_i <- param_list[[1]]
    cov_0_i <- P_i*(1-P_i)
    
    Pr_i <- 1-P_i
    
  } else if (dist == "Pois"){
    lambda_i <- param_list[[1]]
    cov_0_i <- lambda_i
    
    idx_cut_i <- min(which(ppois(0:N,lambda_i)>=1-10^(-8)))
    Pr_i <- ppois(0:(idx_cut_i-2),lambda_i)
    
  }else if (dist == "multinom"){
    P_i <- param_list[[1]]
    num_category <- param_list[[2]]
    mean_i <- P_i%*%c(1:num_category)
    cov_0_i <- (c(1:num_category)-rep(mean_i,num_category))^2%*%P_i
    
    Pr_i <- cumsum(P_i[1:(num_category-1)])
    
  } else if (dist == "negbin"){
    
    P_i <- param_list[[1]]
    r_i <- param_list[[2]]
    cov_0_i <- r_i*(1-P_i)/P_i^2
    
    idx_cut_i <- min(min(which(pnbinom(0:N,r_i,P_i)>=1-10^(-8))),N)
    Pr_i <- pnbinom(0:(idx_cut_i-2),r_i,P_i)
    
  }
  
  # ------------------------------------------------- #
  # jth entry:
  if (dist == "Bern"){
    P_j <- param_list[[2]]
    cov_0_j <- P_j*(1-P_j)
    
    Pr_j <- 1-P_j
    
  } else if (dist == "Pois"){
    lambda_j <- param_list[[2]]
    cov_0_j <- lambda_j
    
    idx_cut_j <- min(which(ppois(0:N,lambda_j)>=1-10^(-8)))
    Pr_j <- ppois(0:(idx_cut_j-2),lambda_j)
    
  }else if (dist == "multinom"){
    P_j <- param_list[[3]]
    num_category <- param_list[[4]]
    mean_j <- P_i%*%c(1:num_category)
    cov_0_j <- (c(1:num_category)-rep(mean_j,num_category))^2%*%P_j
    
    Pr_j <- cumsum(P_j[1:(num_category-1)])
    
  } else if (dist == "negbin"){
    
    P_j <- param_list[[3]]
    r_j <- param_list[[4]]
    cov_0_j <- r_j*(1-P_j)/P_j^2
    
    idx_cut_j <- min(min(which(pnbinom(0:N,r_j,P_j)>=1-10^(-8))),N)
    Pr_j <- pnbinom(0:(idx_cut_j-2),r_j,P_j)
    
  }
  
  # ------------------------------------------------- #
  # Hermite polynomials
  eps <- 1e-8
  if (sum(Pr_i==0)!=0){
    Pr_i[which(Pr_i==0)] <- eps
  }
  if (sum(Pr_i==1)!=0){
    Pr_i[which(Pr_i==1)] <- 1-eps
  }
  if (sum(Pr_j==0)!=0){
    Pr_j[which(Pr_j==0)] <- eps
  }
  if (sum(Pr_j==1)!=0){
    Pr_j[which(Pr_j==1)] <- 1-eps
  }
  
  
  
  Hk_i <- lapply(X=1:K,Y=Pr_i,FUN=function(X,Y) {Her <- as.function(Polys[[X]]); Her(qnorm(Pr_i,0,1))} )
  Hk_j <- lapply(X=1:K,Y=Pr_j,FUN=function(X,Y) {Her <- as.function(Polys[[X]]); Her(qnorm(Pr_j,0,1))} )
  
  # Hermite coefficients
  g_i <- unlist(lapply(X=1:K,FUN = function(X) exp( -qnorm(Pr_i,0,1)^2/2 ) %*% Hk_i[[X]]/ (sqrt(2*pi)*factorial(X))))
  g_j <- unlist(lapply(X=1:K,FUN = function(X) exp( -qnorm(Pr_j,0,1)^2/2 ) %*% Hk_j[[X]]/ (sqrt(2*pi)*factorial(X))))
  
  # Polynomial of link function
  l_ij <- c(factorial(1:K))*c(g_i)*c(g_j)/rep(sqrt(cov_0_i*cov_0_j),K)
  
  return(l_ij)
}


#################################################################################
# Function for adjusting computed link functions
#################################################################################
Link_adjustment <- function(u,coef,knot,dist,param){
  
  if (sum(diff(knot) < 0) > 0){
    
    u_L1 <- which(diff(knot)<0 & knot[-1]>=0)
    u_Lm1 <- which(diff(knot)<0 & knot[-1]<0)
    
    U <- runif(10000,0,1)
    if (dist == "Bern"){
      F1 <- cor(1*(U<param[[1]]),1*(U<param[[2]]))
      Fm1 <- cor(1*(U<param[[1]]),1*((1-U)<param[[2]]))
      
    }else if (dist == "multinom"){
      vec_Fi <- vector("numeric",length(U))
      for (i in 1:(param[[2]])){
        vec_Fi <- vec_Fi + 1*(U <= cumsum(param[[1]])[i])
      }
      vec_Fj <- vector("numeric",length(U))
      for (j in 1:(param[[4]])){
        vec_Fj <- vec_Fj + 1*(U <= cumsum(param[[3]])[j])
      }
      F1 <- cor(vec_Fi,vec_Fj)
      
      vec_Fmi <- vector("numeric",length(U))
      for (i in 1:(param[[2]])){
        vec_Fmi <- vec_Fmi + 1*(U <= cumsum(param[[1]])[i])
      }
      vec_Fmj <- vector("numeric",length(U))
      for (j in 1:(param[[4]])){
        vec_Fmj <- vec_Fmj + 1*((1-U) <= cumsum(param[[3]])[j])
      }
      Fm1 <- cor(vec_Fmi,vec_Fmj)
      
    }else if (dist == "Pois"){
      F1 <- cor(qpois(U,param[[1]]),qpois(U,param[[2]]))
      Fm1 <- cor(qpois(U,param[[1]]),qpois((1-U),param[[2]]))
      
    }else if (dist == "negbin"){
      F1 <- cor(qnbinom(U,param[[2]],param[[1]]),
                qnbinom(U,param[[4]],param[[3]]))
      Fm1 <- cor(qnbinom(U,param[[2]],param[[1]]),
                 qnbinom((1-U),param[[4]],param[[3]]))
      
    }
    
    pow <- 1:length(coef)
    if (length(u_Lm1)!=0){
      u_Lm1 <- c(1,u_Lm1+1)
      for (i in 1:length(u_Lm1)){
        knot[u_Lm1[i]] <- c(coef,-(Fm1 - coef %*% (-1)^pow))%*%(u[u_Lm1[i]]^c(pow,(length(coef)+1)))
      }
    }
    
    if (length(u_L1)!=0){
      u_L1 <- u_L1+1
      for (i in 1:length(u_L1)){
        knot[u_L1[i]] <- c(coef,(F1 - coef %*% (1)^pow))%*%(u[u_L1[i]]^c(pow,(length(coef)+1)))
      }
    }
    
    return(knot)
  }else{
    return(knot)
  }
}




