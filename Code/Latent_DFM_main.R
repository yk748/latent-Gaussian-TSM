#################################################################################
# Main functions for Simulation section
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

# define directories
dir_est <- paste0(getwd(),"/sim_est")
dir_r <- paste0(getwd(),"/sim_r")
dir_forecast <- paste0(getwd(),"/sim_forecast")

#################################################################################
# 1. Estimation
#################################################################################
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:2){ # T

        #---------------------------------------------------#
        # Set up:
        p <- 1
        r <- r_list[j]

        # H: This will be used in forecasting
        H <- 12

        d <- d_list[k]
        TT <- TT_list[l]

        rho <- 0.9
        c <- 0.3

        model_list <- list()
        model_list$Psi <- array(diag(rho,r),dim=c(r,r,p))
        model_list$Cov_eta <- diag(1-rho^2,r)
        model_list$Lambda <- matrix(rnorm(c(d*r),0,1),ncol=r)
        unif <- runif(d,c,1-c)
        model_list$Cov_eps <- diag(unif/(1-unif)*rowSums(model_list$Lambda^2),d)

        dist <- dist_list[i]
        dist_opt <- NULL
        if (i == 1){ #"Bern"
          model_list$Param <- c(rep(0.2,d/3),rep(0.4,d/3),rep(0.7,d/3))

        }else if (i == 2){ #"multinom"
          dist_opt <- 5
          model_list$Param <- c(lapply(X=1:(d/3),FUN=function(X)rep(0.2,dist_opt)),
                                lapply(X=1:(d/3),FUN=function(X)c(0,0.25,0.5,0.25,0)),
                                lapply(X=1:(d/3),FUN=function(X)c(0.45,0,0.1,0,0.45)) )

        }else if (i == 3){ #"Pois"
          model_list$Param <- c(rep(0.1,d/3),rep(1,d/3),rep(10,d/3))

        }else if (i == 4){ #"negbin"
          model_list$Param <- c(lapply(X=1:(d/3),FUN=function(X)c(0.2,3)),
                                lapply(X=1:(d/3),FUN=function(X)c(0.4,3)),
                                lapply(X=1:(d/3),FUN=function(X)c(0.7,3)))

        }

        #---------------------------------------------------#
        # Main loop
        sim_model <- list()
        sim_X <- list()
        sim_Y <- list()
        sim_Z <- list()
        sim_cov_X <- list()
        sim_link <- list()
        sim_cov_Z <- list()
        sim_est <- list()
        sim_time <- list()
        N_sim <- 100
        cat("#----------------------------------# \n")
        cat("Setting ",i,j,k,l," begins \n")
        for (iter in 1:N_sim){

          cat(iter,"th iteration begins \n")
          #---------------------------------------------------#
          # Generate data:
          # Make sure that none of X_{i,t} = 0 for all t
          again <- TRUE
          while (again){
            DGP <- Latent_DFM_Model(p,r,d,(TT+H),dist,dist_opt,model_list,identy_opt=1)

            if (i == 1){
              if (sum(rowSums(DGP$X_t) == 0) == 0 & sum(rowSums(DGP$X_t) == TT) == 0){
                again <- FALSE
              }
            }else if (i == 2){
              possible_num <- c(rep(5,d/3),rep(3,d/3),rep(3,d/3))
              if (sum(mapply(x=1:d,function(x)
                sum(seq(1,5) %in% unique(DGP$X_t[x,])))!=possible_num)==0){
                again <- FALSE
              }
            }else { # i == 3 | i == 4
              if (sum(rowSums(DGP$X_t) == 0) == 0){
                again <- FALSE
              }
            }
          }
          if (iter == 1){
            sim_model$Lambda <- DGP$Lambda
            sim_model$Cov_eps <- DGP$Cov_eps
            sim_model$Psi <- DGP$Psi
            sim_model$Cov_eta <- DGP$Cov_eta
            sim_model$threshold <- DGP$threshold
          }
          sim_X[[iter]] <- DGP$X_t
          sim_Y[[iter]] <- DGP$Y_t
          sim_Z[[iter]] <- DGP$Z_t

          #---------------------------------------------------#
          # Estimation:
          start.time <- proc.time()

          X_t <- DGP$X_t[,-c(TT+1:(TT+H))]
          Cov_X <- Latent_Gauss_Cov(d,TT,p,X_t)
          sim_cov_X[[iter]] <- Cov_X
          Link_param <- Latent_Gauss_Link(d,TT,dist,X_t,dist_opt)
          sim_link[[iter]] <- Link_param
          Cov_Z <- Latent_Gauss_invLink(d,p,Cov_X,Link_param$Link)
          sim_cov_Z[[iter]] <- Cov_Z
          Estim <- Latent_DFM_Estim(r,p,Cov_Z,identy_opt=1,shift=TRUE)
          sim_est[[iter]] <- Estim

          end.time <- proc.time()
          sim_time[[iter]] <- end.time - start.time
          cat("Time taken is",sim_time[[iter]][1],"\n")
        }
        cat("Setting ",i,j,k,l," ends \n")
        cat("#----------------------------------# \n")

        fn <- paste0("Estim","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path <- dir_est
        save(sim_model,sim_X,sim_Y,sim_Z,sim_cov_X,sim_link,sim_cov_Z,sim_est,sim_time,
             file = file.path(path,fn))
      }
    }
  }
}

#################################################################################
# 2. r selection
#################################################################################
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:2){ # T

        #---------------------------------------------------#
        # load file:
        cat("Setting ",i,j,k,l," begins \n")
        load_fn <- paste0("Estim","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_est <- dir_est
        load(paste0(path_est,"/",load_fn))

        #---------------------------------------------------#
        # Set up:
        p <- 1

        H <- 12
        d <- d_list[k]
        TT <- TT_list[l]

        sim_trad_r <- list()
        sim_cv_r <- list()
        N_sim <- 100
        for (iter in 1:N_sim){

          cat(iter,"th iteration begins \n")
          #---------------------------------------------------#
          # estimating r:
          X_t <- sim_X[[iter]][,-c(TT+1:(TT+H))]
          trad_r <- Latent_DFM_r_Selection_trad(d,TT,p,sim_cov_Z[[iter]],identy_opt=1)
          sim_trad_r[[iter]]<- trad_r
          cv_r <- Latent_DFM_r_Selection_BCV(d,TT,p,X_t,sim_link[[iter]]$Link,identy_opt=1,fold=4)
          sim_cv_r[[iter]] <- cv_r
        }
        cat("Setting ",i,j,k,l," ends \n")
        cat("#----------------------------------# \n")

        save_fn <- paste0("r_est","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_r <- dir_r
        save(sim_trad_r,sim_cv_r,file = file.path(path_r,save_fn))
      }
    }
  }
}


#################################################################################
# 3. Forecasting
#################################################################################
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:1){ # T
        
        # #---------------------------------------------------#
        # load file:
        cat("Setting ",i,j,k,l," begins \n")
        load_fn <- paste0("Estim","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_est <- dir_est
        load(paste0(path_est,"/",load_fn))
        
        #---------------------------------------------------#
        # Main loop
        sim_Forecast <- list()
        sim_time_f <- list()
        N_sim <- 100
        cat("#----------------------------------# \n")
        cat("Setting ",i,j,k,l," begins \n")
        for (iter in 1:N_sim){
          
          cat(iter,"th iteration begins \n")
          #---------------------------------------------------#
          # Forecasting:
          start.time <- proc.time()
          
          X_t <- DGP$X_t[,-c(TT+1:(TT+H))]
          # X_t <- sim_X[[iter]][,-c(TT+1:(TT+H))]
          if (i == 1){
            dist <-"Bern"
          }else if (i == 2){
            dist <- "multinom"
          }else if (i == 3){
            dist <- "Pois"
          }else if (i == 4){
            dist <- "negbin"
          }
          Forecast <- Latent_DFM_Forecast(sim_model,p,r,d,TT,X_t,H,dist,dist_opt)
          sim_Forecast[[iter]] <- Forecast
          
          end.time <- proc.time()
          sim_time_f[[iter]] <- end.time - start.time
          cat("Time taken is",sim_time_f[[iter]][1],"\n")
        }
        cat("Setting ",i,j,k,l," ends \n")
        cat("#----------------------------------# \n")
        
        save_fn <- paste0("Forecast","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_forecast <- dir_forecast
        save(sim_model,sim_X,sim_Y,sim_Z,sim_Forecast,sim_time_f,
             file = file.path(path_forecast,save_fn))
      }
    }
  }
}


#################################################################################
# 4. Table & illustration
#################################################################################
#----------------------------------#
# Estimation:
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:2){ # T
        
        #----------------------------------#
        # load file:
        cat("Setting ",i,j,k,l," begins \n")
        load_fn <- paste0("Estim","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_est <- dir_est
        load(paste0(path_est,"/",load_fn))
        
        #----------------------------------#
        # set up:
        dist <- dist_list[i]
        r <- r_list[j]
        d <- d_list[k]
        TT <- TT_list[l]
        dist_opt <- NULL
        Param <- list()
        if (i == 1){ #"Bern"
          Param <- c(rep(0.2,d/3),rep(0.4,d/3),rep(0.7,d/3))
          
        }else if (i == 2){ #"multinom"
          dist_opt <- 5
          Param <- c(lapply(X=1:(d/3),FUN=function(X)rep(0.2,dist_opt)),
                     lapply(X=1:(d/3),FUN=function(X)c(0,0.25,0.5,0.25,0)),
                     lapply(X=1:(d/3),FUN=function(X)c(0.45,0,0.1,0,0.45)) )
          
        }else if (i == 3){ #"Pois"
          Param <- c(rep(0.1,d/3),rep(1,d/3),rep(10,d/3))
          
        }else if (i == 4){ #"negbin"
          Param <- c(lapply(X=1:(d/3),FUN=function(X)c(0.2,3)),
                     lapply(X=1:(d/3),FUN=function(X)c(0.4,3)),
                     lapply(X=1:(d/3),FUN=function(X)c(0.7,3)))
          
        }
        
        #----------------------------------#
        # main loop:
        N_sim <- 100
        tmp_param <- tmp_lambda <- tmp_cov_eps <- tmp_psi <- tmp_cov_eta <- vector("numeric",N_sim) 
        for (iter in 1:N_sim){
          if (i == 1){
            tmp_param[iter] <- sqrt(norm(Param - unlist(sim_link[[iter]]$Param),"2")^2/d)
          }else if (i == 2){
            tmp_param[iter] <- sqrt(norm(unlist(Param) - unlist(sim_link[[iter]]$Param),"2")^2/d)
          }else if (i == 3){
            tmp_param[iter] <- sqrt(norm(Param - unlist(sim_link[[iter]]$Param),"2")^2/d)
          }else if (i == 4){
            tmp_param[iter] <- sqrt(norm(mapply(x=1:d,function(x)Param[[x]][1]) 
                                    - mapply(x=1:d,function(x)sim_link[[iter]]$Param[[x]][1]),"2")^2/d)
          }
          
          tmp_lambda[iter] <- sqrt(norm(sim_model$Lambda - sim_est[[iter]]$Lambda,"F")^2/d)
          tmp_cov_eps[iter] <- sqrt(norm(diag(sim_model$Cov_eps) - sim_est[[iter]]$Cov_eps,"2")^2/d)
          tmp_psi[iter] <- sqrt(norm(matrix(sim_model$Psi,r,r) - matrix(sim_est[[iter]]$Psi,r,r),"F")^2/r)
          tmp_cov_eta[iter] <- sqrt(norm(sim_model$Cov_eta - sim_est[[iter]]$Cov_eta,"F")^2/r)
        }
        
        #----------------------------------#
        # reporting:
        cat("#----------------------------------#","\n")
        cat("Currently, distribution is",dist,"\n")
        cat("the number of factor is",r,"\n")
        cat("the dimension is",d,"\n")
        cat("the sample length is",TT,"\n")
        cat("Parameters are means:",round(mean(tmp_param),4),
            "sds are","(",round(sd(tmp_param),4),")","\n")
        cat("Loadings are means:",round(mean(tmp_lambda),4),
            "sds are","(",round(sd(tmp_lambda),4),")","\n")
        cat("Cov of factor is mean:",round(mean(tmp_cov_eps),4),
            "sd is","(",round(sd(tmp_cov_eps),4),")","\n")
        cat("Transition mat is mean:",round(mean(tmp_psi),4),
            "sd is","(",round(sd(tmp_psi),4),")","\n")
        cat("Cov of VAR is mean:",round(mean(tmp_cov_eta),4),
            "sd is","(",round(sd(tmp_cov_eta),4),")","\n")
        cat("#----------------------------------#","\n")

      }
    }
  }
}


#----------------------------------#
# rank selection:
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

df_total <- data.frame()
for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:2){ # T
        
        #----------------------------------#
        # load file:
        cat("Setting ",i,j,k,l," begins \n")
        load_fn <- paste0("r_est","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_est <- dir_r
        load(paste0(path_est,"/",load_fn))
        
        dist <- dist_list[i]
        r <- r_list[j]
        d <- d_list[k]
        TT <- TT_list[l]
        
        N_sim <- 100
        IC_1 <- unlist(mapply(x=1:N_sim,function(x)sim_trad_r[[x]]$r_IC1))
        IC_2 <- unlist(mapply(x=1:N_sim,function(x)sim_trad_r[[x]]$r_IC2))
        IC_3 <- unlist(mapply(x=1:N_sim,function(x)sim_trad_r[[x]]$r_IC3))
        ED <- unlist(mapply(x=1:N_sim,function(x)sim_trad_r[[x]]$r_ED))
        length_ED <- length(ED)
        if ( length_ED != N_sim){
          ED <- c(ED,rep(0,N_sim-length_ED))
        }
        Fac <- unlist(mapply(x=1:N_sim,function(x)sim_cv_r[[x]]$r_Fac))
        PC <- unlist(mapply(x=1:N_sim,function(x)sim_cv_r[[x]]$r_PC))
        
        df_temp <- data.frame(
          
          dist = rep(dist,N_sim*6),
          num_fac = rep(r,N_sim*6),
          dim = rep(d,N_sim*6),
          length = rep(TT,N_sim*6),
          r_hat = c(IC_1,IC_2,IC_3,ED,Fac,PC),
          methods = c(rep("IC1",N_sim),
                      rep("IC2",N_sim),
                      rep("IC3",N_sim),
                      rep("ED",N_sim),
                      rep("Fac",N_sim),
                      rep("PC",N_sim))
        )
        df_total <- rbind(df_total,df_temp)
        
      }
    }
  }
}

df_total$dim <-factor(df_total$dim,levels=unique(df_total$dim),
                      labels=c("d=15","d=30","d=60","d=90"))
df_total$length <- factor(df_total$length,levels=unique(df_total$length),
                          labels=c("T=100","T=200"))
df_total$num_fac <- factor(df_total$num_fac,levels=unique(df_total$num_fac),
                           labels=c("r=2","r=5"))
df_total$methods <- factor(df_total$methods,levels=unique(df_total$methods))
df_total$dist <- factor(df_total$dist,levels=unique(df_total$dist))

dim.labs <- c("d=15","d=30","d=60","d=90")
names(dim.labs) <- c(TeX("$d=15$"),TeX("$d=30$"),TeX("$d=60$"),TeX("$d=90$"))
length.labs <- c("T=100","T=200")
names(length.labs) <- c(TeX("$T=100$"),TeX("$T=200$"))
fac.labs <- c("r=2","r=5")
names(fac.labs) <- c(TeX("$r=2$"),TeX("$r=5$"))


# Plotting:
pl_bern <- ggplot(df_total[df_total$dist=="Bern",],aes(x=r_hat,fill=methods)) +
  geom_histogram(alpha=0.5,bins=11,position='dodge') +
  scale_x_continuous(breaks=c(0:10)) +
  xlab(TeX("$\\hat{r}$")) +
  ylab("Frequency") + 
  scale_fill_manual(name="",values=rev(colorRampPalette(brewer.pal(9,"Set1"))(6)),
                    labels=levels(df_total$methods)) +
  facet_nested(num_fac~length+dim,
               labeller=labeller(length=length.labs,dim=dim.labs,
                                 fac=fac.labs)) + 
  ggtitle("Bernoulli distributions") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        plot.title = element_text(size=10),
        strip.text.x = element_text(size = 8),
        axis.text=element_text(size=7)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

pl_multi <- ggplot(df_total[df_total$dist=="multinom",],aes(x=r_hat,fill=methods)) +
  geom_histogram(alpha=0.5,bins=11,position='dodge') +
  scale_x_continuous(breaks=c(0:10)) +
  xlab(TeX("$\\hat{r}$")) +
  ylab("Frequency") + 
  scale_fill_manual(name="",values=rev(colorRampPalette(brewer.pal(9,"Set1"))(6)),
                    labels=levels(df_total$methods)) +
  facet_nested(num_fac~length+dim,
               labeller=labeller(length=length.labs,dim=dim.labs,
                                 fac=fac.labs)) + 
  ggtitle("Multinomial distributions") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        plot.title = element_text(size=10),
        strip.text.x = element_text(size = 8),
        axis.text=element_text(size=7)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

pl_pois <- ggplot(df_total[df_total$dist=="Pois",],aes(x=r_hat,fill=methods)) +
  geom_histogram(alpha=0.5,bins=11,position='dodge') +
  scale_x_continuous(breaks=c(0:10)) +
  xlab(TeX("$\\hat{r}$")) +
  ylab("Frequency") + 
  scale_fill_manual(name="",values=rev(colorRampPalette(brewer.pal(9,"Set1"))(6)),
                    labels=levels(df_total$methods)) +
  facet_nested(num_fac~length+dim,
               labeller=labeller(length=length.labs,dim=dim.labs,
                                 fac=fac.labs)) + 
  ggtitle("Poisson distributions") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        plot.title = element_text(size=10),
        strip.text.x = element_text(size = 8),
        axis.text=element_text(size=7)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

pl_negbin <- ggplot(df_total[df_total$dist=="negbin",],aes(x=r_hat,fill=methods)) +
  geom_histogram(alpha=0.5,bins=11,position='dodge') +
  scale_x_continuous(breaks=c(0:10)) +
  xlab(TeX("$\\hat{r}$")) +
  ylab("Frequency") + 
  scale_fill_manual(name="",values=rev(colorRampPalette(brewer.pal(9,"Set1"))(6)),
                    labels=levels(df_total$methods)) +
  facet_nested(num_fac~length+dim,
               labeller=labeller(length=length.labs,dim=dim.labs,
                                 fac=fac.labs)) + 
  ggtitle("Negative binomial distributions") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        plot.title = element_text(size=10),
        strip.text.x = element_text(size = 8),
        axis.text=element_text(size=7)) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE))


pdf(width=9,height=12,file = "rank_estim.pdf")
par(mar=c(0,0,0,0))
ggarrange(pl_bern,pl_multi, pl_pois,pl_negbin,
          ncol=1, nrow=4, common.legend = TRUE, legend="bottom")
dev.off()



#----------------------------------#
# forecasting:
dist_list <- c("Bern","multinom","Pois","negbin")
r_list <- c(2,5)
d_list <- c(15,30,60,90)
TT_list <- c(100,200)

for (i in 1:4){ # dist
  for (j in 1:2){ # rank
    for (k in 1:4){ # d
      for (l in 1:1){ # T
        
        #----------------------------------#
        # load file:
        cat("Setting ",i,j,k,l," begins \n")
        load_fn <- paste0("Forecast","_dist",i,"_type",j,"_d",k,"_T",l,".rda")
        path_est <- dir_forecast
        load(paste0(path_est,"/",load_fn))
        
        #----------------------------------#
        # set up:
        dist <- dist_list[i]
        r <- r_list[j]
        d <- d_list[k]
        TT <- TT_list[l]
        
        H1 <- 1
        H2 <- 2
        H3 <- 3
        H6 <- 6
        H12 <- 12
        
        N_sim <- 100
        N <- 100
        Y_H1 <- Y_H2 <- Y_H3 <- Y_H6 <- Y_H12 <- vector("numeric",N_sim) 
        Z_H1 <- Z_H2 <- Z_H3 <- Z_H6 <- Z_H12 <- vector("numeric",N_sim) 
        X_H1 <- X_H2 <- X_H3 <- X_H6 <- X_H12 <- vector("numeric",N_sim) 
        ACC_H1 <- ACC_H2 <- ACC_H3 <- ACC_H6 <- ACC_H12 <- vector("numeric",N_sim)
        L_H1 <- L_H2 <- L_H3 <- L_H6 <- L_H12 <- vector("numeric",N_sim)
        M_H1 <- M_H2 <- M_H3 <- M_H6 <- M_H12 <- vector("numeric",N_sim)
        I_H1 <- I_H2 <- I_H3 <- I_H6 <- I_H12 <- vector("numeric",N_sim)
        for (iter in 1:N_sim){
          Y_H1[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Y[[iter]][,TT+H1] 
                                                                 - sim_Forecast[[iter]]$Y_hat_h[,H1,x],"2")^2))/r)
          Y_H2[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Y[[iter]][,TT+H2] 
                                                                 - sim_Forecast[[iter]]$Y_hat_h[,H2,x],"2")^2))/r)
          Y_H3[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Y[[iter]][,TT+H3] 
                                                                 - sim_Forecast[[iter]]$Y_hat_h[,H3,x],"2")^2))/r)
          Y_H6[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Y[[iter]][,TT+H6] 
                                                                 - sim_Forecast[[iter]]$Y_hat_h[,H6,x],"2")^2))/r)
          Y_H12[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Y[[iter]][,TT+H12] 
                                                                  - sim_Forecast[[iter]]$Y_hat_h[,H12,x],"2")^2))/r)
          
          Z_H1[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Z[[iter]][,TT+H1] 
                                                                 - sim_Forecast[[iter]]$Z_hat_h[,H1,x],"2")^2))/d)
          Z_H2[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Z[[iter]][,TT+H2] 
                                                                 - sim_Forecast[[iter]]$Z_hat_h[,H2,x],"2")^2))/d)
          Z_H3[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Z[[iter]][,TT+H3] 
                                                                 - sim_Forecast[[iter]]$Z_hat_h[,H3,x],"2")^2))/d)
          Z_H6[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Z[[iter]][,TT+H6] 
                                                                 - sim_Forecast[[iter]]$Z_hat_h[,H6,x],"2")^2))/d)
          Z_H12[iter] <- sqrt(mean(mapply(x=1:100,function(x)norm(sim_Z[[iter]][,TT+H12] 
                                                                  - sim_Forecast[[iter]]$Z_hat_h[,H12,x],"2")^2))/d)
          
          X_H1[iter] <- sqrt(norm(sim_X[[iter]][,TT+H1] - sim_Forecast[[iter]]$X_t_h[,H1],"2")/d)
          X_H2[iter] <- sqrt(norm(sim_X[[iter]][,TT+H2] - sim_Forecast[[iter]]$X_t_h[,H2],"2")/d)
          X_H3[iter] <- sqrt(norm(sim_X[[iter]][,TT+H3] - sim_Forecast[[iter]]$X_t_h[,H3],"2")/d)
          X_H6[iter] <- sqrt(norm(sim_X[[iter]][,TT+H6] - sim_Forecast[[iter]]$X_t_h[,H6],"2")/d)
          X_H12[iter] <- sqrt(norm(sim_X[[iter]][,TT+H12] - sim_Forecast[[iter]]$X_t_h[,H12],"2")/d)
          
          ACC_H1[iter] <- sum(sim_X[[iter]][,TT+H1]==sim_Forecast[[iter]]$X_t_h[,H1])/d
          ACC_H2[iter] <- sum(sim_X[[iter]][,TT+H2]==sim_Forecast[[iter]]$X_t_h[,H2])/d
          ACC_H3[iter] <- sum(sim_X[[iter]][,TT+H3]==sim_Forecast[[iter]]$X_t_h[,H3])/d
          ACC_H6[iter] <- sum(sim_X[[iter]][,TT+H6]==sim_Forecast[[iter]]$X_t_h[,H6])/d
          ACC_H12[iter] <- sum(sim_X[[iter]][,TT+H12]==sim_Forecast[[iter]]$X_t_h[,H12])/d
          
          L <- sim_X[[iter]][,TT]
          L_H1[iter] <- sum(sim_X[[iter]][,TT+H1] == L)/d
          L_H2[iter] <- sum(sim_X[[iter]][,TT+H2] == L)/d
          L_H3[iter] <- sum(sim_X[[iter]][,TT+H3] == L)/d
          L_H6[iter] <- sum(sim_X[[iter]][,TT+H6] == L)/d
          L_H12[iter] <- sum(sim_X[[iter]][,TT+H12] == L)/d
          
          M <- mapply(x=1:d,function(x)data.frame(sort(table(sim_X[[iter]][x,1:TT]),decreasing=TRUE))[1,1])
          M_H1[iter] <- sum(sim_X[[iter]][,TT+H1] == M)/d
          M_H2[iter] <- sum(sim_X[[iter]][,TT+H2] == M)/d
          M_H3[iter] <- sum(sim_X[[iter]][,TT+H3] == M)/d
          M_H6[iter] <- sum(sim_X[[iter]][,TT+H6] == M)/d
          M_H12[iter] <- sum(sim_X[[iter]][,TT+H12] == M)/d
          
          if (i == 1){
            Z0 <- 1*(unlist(sim_model$threshold) > 0)
          }else if (i == 2){
            Z0 <- mapply(x=1:d,function(x){sum(sim_model$threshold[[x]] < 0)+1})
          }else{
            Z0 <- mapply(x=1:d,function(x){sum(sim_model$threshold[[x]] < 0)})
          }
          I_H1[iter] <- sum(sim_X[[iter]][,TT+H1] == Z0)/d
          I_H2[iter] <- sum(sim_X[[iter]][,TT+H2] == Z0)/d
          I_H3[iter] <- sum(sim_X[[iter]][,TT+H3] == Z0)/d
          I_H6[iter] <- sum(sim_X[[iter]][,TT+H6] == Z0)/d
          I_H12[iter] <- sum(sim_X[[iter]][,TT+H12] == Z0)/d
        }

        #----------------------------------#
        # reporting:
        cat("#----------------------------------#","\n")
        cat("Currently, distribution is",dist,"\n")
        cat("the number of factor is",r,"\n")
        cat("the dimension is",d,"\n")
        cat("the sample length is",TT,"\n")
        cat("MSFE_Y at H=1 are mean (sd):",round(mean(Y_H1),4),"(",round(sd(Y_H1),4),")","\n")
        cat("MSFE_Y at H=2 are mean (sd):",round(mean(Y_H2),4),"(",round(sd(Y_H2),4),")","\n")
        cat("MSFE_Y at H=3 are mean (sd):",round(mean(Y_H3),4),"(",round(sd(Y_H3),4),")","\n")
        cat("MSFE_Y at H=6 are mean (sd):",round(mean(Y_H6),4),"(",round(sd(Y_H6),4),")","\n")
        cat("MSFE_Y at H=12 are mean (sd):",round(mean(Y_H12),4),"(",round(sd(Y_H12),4),")","\n")
        cat("MSFE_Z at H=1 are mean (sd):",round(mean(Z_H1),4),"(",round(sd(Z_H1),4),")","\n")
        cat("MSFE_Z at H=2 are mean (sd):",round(mean(Z_H2),4),"(",round(sd(Z_H2),4),")","\n")
        cat("MSFE_Z at H=3 are mean (sd):",round(mean(Z_H3),4),"(",round(sd(Z_H3),4),")","\n")
        cat("MSFE_Z at H=6 are mean (sd):",round(mean(Z_H6),4),"(",round(sd(Z_H6),4),")","\n")
        cat("MSFE_Z at H=12 are mean (sd):",round(mean(Z_H12),4),"(",round(sd(Z_H12),4),")","\n")
        cat("MSFE_X at H=1 are mean (sd):",round(mean(X_H1),4),"\n")
        cat("MSFE_X at H=2 are mean (sd):",round(mean(X_H2),4),"\n")
        cat("MSFE_X at H=3 are mean (sd):",round(mean(X_H3),4),"\n")
        cat("MSFE_X at H=6 are mean (sd):",round(mean(X_H6),4),"\n")
        cat("MSFE_X at H=12 are mean (sd):",round(mean(X_H12),4),"\n")
        cat("ACC at H=1 are mean (L,M,I):",round(mean(ACC_H1),4),"(",round(c(mean(L_H1),mean(M_H1),mean(I_H1)),4),")","\n")
        cat("ACC at H=2 are mean (L,M,I):",round(mean(ACC_H2),4),"(",round(c(mean(L_H2),mean(M_H2),mean(I_H2)),4),")","\n")
        cat("ACC at H=3 are mean (L,M,I):",round(mean(ACC_H3),4),"(",round(c(mean(L_H3),mean(M_H3),mean(I_H3)),4),")","\n")
        cat("ACC at H=6 are mean (L,M,I):",round(mean(ACC_H6),4),"(",round(c(mean(L_H6),mean(M_H6),mean(I_H6)),4),")","\n")
        cat("ACC at H=12 are mean (L,M,I):",round(mean(ACC_H12),4),"(",round(c(mean(L_H12),mean(M_H12),mean(I_H12)),4),")","\n")

        cat("#----------------------------------#","\n")
      }
    }
  }
}


