#################################################################################
# Main functions for Application section
#################################################################################
rm(list=ls())

# load packages:
library("polynom")
library("matlib")
library("latex2exp")
library("ggplot2")
library("patchwork")
library("ggpubr")
library("RColorBrewer")
library("reshape2")
library("ggh4x")
library("mvtnorm")
library("MASS")
library("TruncatedNormal")
library("psych")
library("GPArotation")
library("Matrix")

# load source files:
load(file = "Polys.rda")
source("Latent_Gauss_Cov.R")
source("Latent_Gauss_Link.R")
source("Latent_Gauss_invLink.R")
source("Latent_DFM_Model.R")
source("Latent_DFM_Estim.R")
source("Latent_DFM_r_Selection_trad.R")
source("Latent_DFM_r_Selection_BCV.R")
source("Latent_DFM_Forecast.R")


#################################################################################
# Estimation:
#################################################################################
# set up:
load("Data.rda")

d <- 30
TT <- 85
dist <- "multinom"
dist_opt <- 5
p <- 1
H <- 5

X_t <- data[,-c(TT+1:(TT+H))]
#---------------------------------------------------#
# Link function:
Link_param <- Latent_Gauss_Link(d,TT,dist,X_t,dist_opt)

#---------------------------------------------------#
# Estimation:
Cov_X <- Latent_Gauss_Cov(d,TT,p,X_t)
Cov_Z <- Latent_Gauss_invLink(d,p,Cov_X,Link_param$Link)
# Cov_Z <- Latent_Gauss_invLink(d,p,Cov_X,Link_param$Link,dist="multinom",Link_param$Param)
Estim <- Latent_DFM_Estim(r=5,p,Cov_Z,identy_opt=2,shift=TRUE)


# ------------------------------------------------- #
# Set up from models
r <- 5
Lambda <- Estim$Lambda
Cov_eps = Estim$Cov_eps
Psi <- companion_form_phi(Estim$Psi,r,p)
Cov_eta <- Estim$Cov_eta
Param <- Link_param$Param

# ------------------------------------------------- #
# Set up for Kalman recursions
window_sample <- TT
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
threshold <- NULL
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
  # R_hat[,,t] <- Lambda %*% Q_hat[,,t] %*% t(Lambda) + Cov_eps
  R_hat[,,t] <- Lambda %*% Q_hat[,,t] %*% t(Lambda) + diag(abs(Estim$Cov_eps),d)
  
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
# Compute PIT
# # 1st variable.


P_t_y <- array(NA,dim=c(d,TT,5))
for (i in 1:d){
  aug_thres <- c(-Inf,threshold[[i]],Inf)
  for (t in 1:TT){
    for (l in 1:5){
      P_t_y[i,t,l] <- sum(w_t[t,] * mapply(k=1:100,function(k)pnorm((aug_thres[(l+1)] - Z_hat[i,t,k])/R_hat[i,i,t]) - pnorm((aug_thres[l] - Z_hat[i,t,k])/R_hat[i,i,t])))
    }
    P_t_y[i,t,] <- cumsum(P_t_y[i,t,])
  }
}


HH <- 11
F_t_y <- array(NA,dim=c(d,TT,HH))
F_bar <- array(NA,dim=c(d,HH))
for (i in 1:d){
  for (hh in 1:HH){
    u <- (hh-1)/HH
    
    for (t in 1:TT){
      if (X_t[i,t] == 1){
        lower <- 0
        upper <- P_t_y[i,t,1]
      }else{
        lower <- P_t_y[i,t,(X_t[i,t]-1)]
        upper <- P_t_y[i,t,X_t[i,t]]
      }
      
      if (upper == lower){
        F_t_y[i,t,] <- 0
      }else{
        if (u<=lower){
          F_t_y[i,t,hh] <- 0
        }else if (u>=upper){
          F_t_y[i,t,hh] <- 1
        }else{
          F_t_y[i,t,hh] <- (u - lower)/(upper - lower)
        }
      }
    }
    F_bar[i,hh] <- sum(F_t_y[i,,hh])/(TT+1)
  }
}


heights <- mapply(i=1:d,function(i)diff(F_bar[i,])/sum(diff(F_bar[i,])))
df_F_bar <- data.frame(t(heights))
colnames(df_F_bar) <- (c(1:(HH-1)) - 1)/(HH-1)
rownames(df_F_bar) <- rownames(data)
melt_df_F_bar <- melt(t(df_F_bar))

Mycol <- c(brewer.pal(n=7,name="Blues")[-1],
           brewer.pal(n=7,name="Greens")[-1],
           brewer.pal(n=7,name="Oranges")[-1],
           brewer.pal(n=7,name="Purples")[-1],
           brewer.pal(n=7,name="Greys")[-1])


pdf(width=12,height=8,file="PIT.pdf")
ggplot(melt_df_F_bar, aes(x = Var1, y=value, fill=factor(Var2))) +
  geom_col(color = "black", alpha = 0.6) +
  geom_hline(yintercept = 0.1, linetype="dashed", color="black",linewidth=0.3) + 
  scale_fill_manual(values=Mycol) + 
  scale_x_continuous(label=c(0,0.5,1),breaks=c(-0.05,0.45,0.95)) +
  facet_wrap(~ Var2, scales = "free_y") +  
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "",
       x = "",
       y = "Relative Frequency")
dev.off()
