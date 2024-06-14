#################################################################################
# Function for illustrating Figure 1
#################################################################################
rm(list=ls())

# load packages:
library("polynom")
library("matlib")
library("latex2exp")
library("ggplot2")
library("patchwork")
library("reshape2")

# load source files:
load(file = "Polys.rda")
source("Latent_Gauss_Cov.R")
source("Latent_Gauss_Link.R")
source("Latent_Gauss_invLink.R")

#################################################################################
# 0. Illustration - Link function and inverse link function
#################################################################################
u <- c(seq(-1,-0.72,0.01),seq(-0.7,0,0.2),seq(0.1,0.7,0.2),seq(0.72,1,0.01)) # Pre-set grid

#---------------------------------------------------#
# Bernoulli distribution:
# Set parameters:
param_case1 <- list()
param_case1[[1]] <- 0.4
param_case1[[2]] <- 0.4
l_ij_case1 <- l_ij_func(param_case1,dist="Bern",K=100)

param_case2 <- list()
param_case2[[1]] <- 0.2
param_case2[[2]] <- 0.7
l_ij_case2 <- l_ij_func(param_case2,dist="Bern",K=100)

param_case3 <- list()
param_case3[[1]] <- 0.2
param_case3[[2]] <- 0.4
l_ij_case3 <- l_ij_func(param_case3,dist="Bern",K=100)

param_case4 <- list()
param_case4[[1]] <- 0.4
param_case4[[2]] <- 0.7
l_ij_case4 <- l_ij_func(param_case4,dist="Bern",K=100)

# Compute link and inverse link functions:
L_case1 <- L_case2 <- L_case3 <- L_case4 <- vector("numeric",length=length(u))
for (i in 1:length(u)){
  pow <- 1:length(l_ij_case1)
  L_case1[i] <- l_ij_case1[1:length(l_ij_case1)]%*%(u[i]^pow)[1:length(l_ij_case1)]
  L_case2[i] <- l_ij_case2[1:length(l_ij_case2)]%*%(u[i]^pow)[1:length(l_ij_case2)]
  L_case3[i] <- l_ij_case3[1:length(l_ij_case3)]%*%(u[i]^pow)[1:length(l_ij_case3)]
  L_case4[i] <- l_ij_case4[1:length(l_ij_case4)]%*%(u[i]^pow)[1:length(l_ij_case4)]
}

v <- seq(-1,1,0.01) # Assumed input Cov_X
L_inv_case1 <- interpolation(l_ij_case1,u,v)
L_inv_case2 <- interpolation(l_ij_case2,u,v)
L_inv_case3 <- interpolation(l_ij_case3,u,v)
L_inv_case4 <- interpolation(l_ij_case4,u,v)

# plotting:
df_L <- data.frame(u=u,case1=L_case1,case2=L_case2,case3=L_case3,case4=L_case4)
df_L_melt <- melt(df_L,id="u")

df_inv_L1 <- data.frame(x=v,case1=L_inv_case1,case2=L_inv_case2,case3=L_inv_case3,case4=L_inv_case4)
df_inv_L1_melt <- melt(df_inv_L1,id="x")
df_inv_L1_melt$type <- "interpolated"
df_inv_L2_melt <- data.frame(x=c(L_case1,L_case2,L_case3,L_case4),
                             variable=c(rep("case1",length(u)),rep("case2",length(u)),
                                        rep("case3",length(u)),rep("case4",length(u))),
                             value=rep(u,4))
df_inv_L2_melt$type <- "exact"
df_inv_L_melt <- rbind(df_inv_L1_melt,df_inv_L2_melt)
df_inv_L_melt$type <-factor(df_inv_L_melt$type,levels=c("interpolated","exact"))

cols <- unname(TeX(c("$p_{i}=0.4,p_{j}=0.4$",
                     "$p_{i}=0.2,p_{j}=0.7$", 
                     "$p_{i}=0.2,p_{j}=0.4$",
                     "$p_{i}=0.4,p_{j}=0.7$")))
lines <- c("Interpolated","Exact")

pdf(width=9.5,height=9.5,file = "link_Bern.pdf")
par(mar=c(0,0,0,0))
ggplot(df_L_melt,aes(y=value, x=u, color=variable)) +
  scale_y_continuous(name=TeX("$L(u)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$u$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("Bern($p_i$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  theme(legend.position=c(0.75,0.2),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  geom_line(cex=0.5)
dev.off()

pdf(width=9.5,height=9.5,file = "inv_link_Bern.pdf")
par(mar=c(0,0,0,0))
ggplot(df_inv_L_melt,aes(y=value, x=x, color=variable, linetype=type)) +
  scale_y_continuous(name=TeX("$L^{-1}(v)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$v$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("Bern($p_i$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  scale_linetype_manual(name="",
                     labels=lines,values=c("dashed","dotted")) +
  theme(legend.position=c(0.75,0.3),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  guides(color = guide_legend(order=1), linetype = guide_legend(order=0)) +
  geom_line(cex=0.5)
dev.off()

#---------------------------------------------------#
# Multinomial distribution:
# Set parameters:
param_case1 <- list()
param_case1[[1]] <- rep(0.2,5)
param_case1[[2]] <- 5
param_case1[[3]] <- rep(0.2,5)
param_case1[[4]] <- 5
l_ij_case1 <- l_ij_func(param_case1,dist="multinom",K=100)

param_case2 <- list()
param_case2[[1]] <- rep(0.2,5)
param_case2[[2]] <- 5
param_case2[[3]] <- c(0.45,0,0.1,0,0.45)
param_case2[[4]] <- 5
l_ij_case2 <- l_ij_func(param_case2,dist="multinom",K=100)

param_case3 <- list()
param_case3[[1]] <- rep(0.2,5)
param_case3[[2]] <- 5
param_case3[[3]] <- c(0,0.25,0.5,0.25,0)
param_case3[[4]] <- 5
l_ij_case3 <- l_ij_func(param_case3,dist="multinom",K=100)

param_case4 <- list()
param_case4[[1]] <- c(0,0.25,0.5,0.25,0)
param_case4[[2]] <- 5
param_case4[[3]] <- c(0.45,0,0.1,0,0.45)
param_case4[[4]] <- 5
l_ij_case4 <- l_ij_func(param_case4,dist="multinom",K=100)

# Compute link and inverse link functions:
L_case1 <- L_case2 <- L_case3 <- L_case4 <- vector("numeric",length=length(u))
for (i in 1:length(u)){
  pow <- 1:length(l_ij_case1)
  L_case1[i] <- l_ij_case1[1:length(l_ij_case1)]%*%(u[i]^pow)[1:length(l_ij_case1)]
  L_case2[i] <- l_ij_case2[1:length(l_ij_case2)]%*%(u[i]^pow)[1:length(l_ij_case2)]
  L_case3[i] <- l_ij_case3[1:length(l_ij_case3)]%*%(u[i]^pow)[1:length(l_ij_case3)]
  L_case4[i] <- l_ij_case4[1:length(l_ij_case4)]%*%(u[i]^pow)[1:length(l_ij_case4)]
}

v <- seq(-1,1,0.01) # Assumed input Cov_X
L_inv_case1 <- interpolation(l_ij_case1,u,v)
L_inv_case2 <- interpolation(l_ij_case2,u,v)
L_inv_case3 <- interpolation(l_ij_case3,u,v)
L_inv_case4 <- interpolation(l_ij_case4,u,v)

# plotting:
df_L <- data.frame(u=u,case1=L_case1,case2=L_case2,case3=L_case3,case4=L_case4)
df_L_melt <- melt(df_L,id="u")

df_inv_L1 <- data.frame(x=v,case1=L_inv_case1,case2=L_inv_case2,case3=L_inv_case3,case4=L_inv_case4)
df_inv_L1_melt <- melt(df_inv_L1,id="x")
df_inv_L1_melt$type <- "interpolated"
df_inv_L2_melt <- data.frame(x=c(L_case1,L_case2,L_case3,L_case4),
                             variable=c(rep("case1",length(u)),rep("case2",length(u)),
                                        rep("case3",length(u)),rep("case4",length(u))),
                             value=rep(u,4))
df_inv_L2_melt$type <- "exact"
df_inv_L_melt <- rbind(df_inv_L1_melt,df_inv_L2_melt)
df_inv_L_melt$type <-factor(df_inv_L_melt$type,levels=c("interpolated","exact"))

cols <- unname(TeX(c("$p_{i}=\t{uniform},p_{j}=\t{uniform}$",
                     "$p_{i}=\t{uniform},p_{j}=\t{trimodal}$", 
                     "$p_{i}=\t{unimodal},p_{j}=\t{uniform}$",
                     "$p_{i}=\t{unimodal},p_{j}=\t{trimodal}$")))
lines <- c("Interpolated","Exact")

pdf(width=9.5,height=9.5,file = "link_multinom.pdf")
par(mar=c(0,0,0,0))
ggplot(df_L_melt,aes(y=value, x=u, color=variable)) +
  scale_y_continuous(name=TeX("$L(u)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$u$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("multinom($p_{i}=(p_{i1},...,p_{i5})$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  theme(legend.position=c(0.75,0.2),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  geom_line(cex=0.5)
dev.off()

pdf(width=9.5,height=9.5,file = "inv_link_multinom.pdf")
par(mar=c(0,0,0,0))
ggplot(df_inv_L_melt,aes(y=value, x=x, color=variable, linetype=type)) +
  scale_y_continuous(name=TeX("$L^{-1}(v)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$v$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("multinom($p_{i}=(p_{i1},...,p_{i5})$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  scale_linetype_manual(name="",
                        labels=lines,values=c("dashed","dotted")) +
  theme(legend.position=c(0.78,0.28),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +          guides(color = guide_legend(order=1), linetype = guide_legend(order=0)) +
  geom_line(cex=0.5)
dev.off()

#---------------------------------------------------#
# Poisson distribution:
# Set parameters:
param_case1 <- list()
param_case1[[1]] <- 1
param_case1[[2]] <- 1
l_ij_case1 <- l_ij_func(param_case1,dist="Pois",K=100)

param_case2 <- list()
param_case2[[1]] <- 0.1
param_case2[[2]] <- 10
l_ij_case2 <- l_ij_func(param_case2,dist="Pois",K=100)

param_case3 <- list()
param_case3[[1]] <- 0.1
param_case3[[2]] <- 1
l_ij_case3 <- l_ij_func(param_case3,dist="Pois",K=100)

param_case4 <- list()
param_case4[[1]] <- 1
param_case4[[2]] <- 10
l_ij_case4 <- l_ij_func(param_case4,dist="Pois",K=100)

# Compute link and inverse link functions:
L_case1 <- L_case2 <- L_case3 <- L_case4 <- vector("numeric",length=length(u))
for (i in 1:length(u)){
  pow <- 1:length(l_ij_case1)
  L_case1[i] <- l_ij_case1[1:length(l_ij_case1)]%*%(u[i]^pow)[1:length(l_ij_case1)]
  L_case2[i] <- l_ij_case2[1:length(l_ij_case2)]%*%(u[i]^pow)[1:length(l_ij_case2)]
  L_case3[i] <- l_ij_case3[1:length(l_ij_case3)]%*%(u[i]^pow)[1:length(l_ij_case3)]
  L_case4[i] <- l_ij_case4[1:length(l_ij_case4)]%*%(u[i]^pow)[1:length(l_ij_case4)]
}

v <- seq(-1,1,0.01) # Assumed input Cov_X
L_inv_case1 <- interpolation(l_ij_case1,u,v)
L_inv_case2 <- interpolation(l_ij_case2,u,v)
L_inv_case3 <- interpolation(l_ij_case3,u,v)
L_inv_case4 <- interpolation(l_ij_case4,u,v)


# plotting:
df_L <- data.frame(u=u,case1=L_case1,case2=L_case2,case3=L_case3,case4=L_case4)
df_L_melt <- melt(df_L,id="u")

df_inv_L1 <- data.frame(x=v,case1=L_inv_case1,case2=L_inv_case2,case3=L_inv_case3,case4=L_inv_case4)
df_inv_L1_melt <- melt(df_inv_L1,id="x")
df_inv_L1_melt$type <- "interpolated"
df_inv_L2_melt <- data.frame(x=c(L_case1,L_case2,L_case3,L_case4),
                             variable=c(rep("case1",length(u)),rep("case2",length(u)),
                                        rep("case3",length(u)),rep("case4",length(u))),
                             value=rep(u,4))
df_inv_L2_melt$type <- "exact"
df_inv_L_melt <- rbind(df_inv_L1_melt,df_inv_L2_melt)
df_inv_L_melt$type <-factor(df_inv_L_melt$type,levels=c("interpolated","exact"))

cols <- unname(TeX(c("$\\theta_{i}=1,\\theta_{j}=1$",
                     "$\\theta_{i}=0.1,\\theta_{j}=10$", 
                     "$\\theta_{i}=0.1,\\theta_{j}=1$",
                     "$\\theta_{i}=1,\\theta_{j}=10$")))
lines <- c("Interpolated","Exact")

pdf(width=9.5,height=9.5,file = "link_Pois.pdf")
par(mar=c(0,0,0,0))
ggplot(df_L_melt,aes(y=value, x=u, color=variable)) +
  scale_y_continuous(name=TeX("$L(u)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$u$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("Pois($\\theta_{i}$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  theme(legend.position=c(0.75,0.2),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  geom_line(cex=0.5)
dev.off()

pdf(width=9.5,height=9.5,file = "inv_link_Pois.pdf")
par(mar=c(0,0,0,0))
ggplot(df_inv_L_melt,aes(y=value, x=x, color=variable, linetype=type)) +
  scale_y_continuous(name=TeX("$L^{-1}(v)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$v$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("Pois($\\theta_{i}$)"),
                     labels=cols,values=c("black","blue","red","green")) +
  scale_linetype_manual(name="",
                        labels=lines,values=c("dashed","dotted")) +
  theme(legend.position=c(0.75,0.3),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  guides(color = guide_legend(order=1), linetype = guide_legend(order=0)) +
  geom_line(cex=0.5)
dev.off()

#---------------------------------------------------#
# Negative binomial distribution:
# Set parameters:
param_case1 <- list()
param_case1[[1]] <- 0.4
param_case1[[2]] <- 3
param_case1[[3]] <- 0.4
param_case1[[4]] <- 3
l_ij_case1 <- l_ij_func(param_case1,dist="negbin",K=100)

param_case2 <- list()
param_case2[[1]] <- 0.2
param_case2[[2]] <- 3
param_case2[[3]] <- 0.7
param_case2[[4]] <- 3
l_ij_case2 <- l_ij_func(param_case2,dist="negbin",K=100)

param_case3 <- list()
param_case3[[1]] <- 0.2
param_case3[[2]] <- 3
param_case3[[3]] <- 0.4
param_case3[[4]] <- 3
l_ij_case3 <- l_ij_func(param_case3,dist="negbin",K=100)

param_case4 <- list()
param_case4[[1]] <- 0.4
param_case4[[2]] <- 3
param_case4[[3]] <- 0.7
param_case4[[4]] <- 3
l_ij_case4 <- l_ij_func(param_case4,dist="negbin",K=100)

# Compute link and inverse link functions:
L_case1 <- L_case2 <- L_case3 <- L_case4 <- vector("numeric",length=length(u))
for (i in 1:length(u)){
  pow <- 1:length(l_ij_case1)
  L_case1[i] <- l_ij_case1[1:length(l_ij_case1)]%*%(u[i]^pow)[1:length(l_ij_case1)]
  L_case2[i] <- l_ij_case2[1:length(l_ij_case2)]%*%(u[i]^pow)[1:length(l_ij_case2)]
  L_case3[i] <- l_ij_case3[1:length(l_ij_case3)]%*%(u[i]^pow)[1:length(l_ij_case3)]
  L_case4[i] <- l_ij_case4[1:length(l_ij_case4)]%*%(u[i]^pow)[1:length(l_ij_case4)]
}

v <- seq(-1,1,0.01) # Assumed input Cov_X
L_inv_case1 <- interpolation(l_ij_case1,u,v)
L_inv_case2 <- interpolation(l_ij_case2,u,v)
L_inv_case3 <- interpolation(l_ij_case3,u,v)
L_inv_case4 <- interpolation(l_ij_case4,u,v)


# plotting:
df_L <- data.frame(u=u,case1=L_case1,case2=L_case2,case3=L_case3,case4=L_case4)
df_L_melt <- melt(df_L,id="u")

df_inv_L1 <- data.frame(x=v,case1=L_inv_case1,case2=L_inv_case2,case3=L_inv_case3,case4=L_inv_case4)
df_inv_L1_melt <- melt(df_inv_L1,id="x")
df_inv_L1_melt$type <- "interpolated"
df_inv_L2_melt <- data.frame(x=c(L_case1,L_case2,L_case3,L_case4),
                             variable=c(rep("case1",length(u)),rep("case2",length(u)),
                                        rep("case3",length(u)),rep("case4",length(u))),
                             value=rep(u,4))
df_inv_L2_melt$type <- "exact"
df_inv_L_melt <- rbind(df_inv_L1_melt,df_inv_L2_melt)
df_inv_L_melt$type <-factor(df_inv_L_melt$type,levels=c("interpolated","exact"))

cols <- unname(TeX(c("$p_{i}=0.4,p_{j}=0.4$",
                     "$p_{i}=0.2,p_{j}=0.7$", 
                     "$p_{i}=0.2,p_{j}=0.4$",
                     "$p_{i}=0.4,p_{j}=0.7$")))
lines <- c("Interpolated","Exact")

pdf(width=9.5,height=9.5,file = "link_negbin.pdf")
par(mar=c(0,0,0,0))
ggplot(df_L_melt,aes(y=value, x=u, color=variable)) +
  scale_y_continuous(name=TeX("$L(u)$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$u$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("NB$(p_{i},r=3)$"),
                     labels=cols,values=c("black","blue","red","green")) +
  theme(legend.position=c(0.75,0.2),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  geom_line(cex=0.5)
dev.off()

pdf(width=9.5,height=9.5,file = "inv_link_negbin.pdf")
par(mar=c(0,0,0,0))
ggplot(df_inv_L_melt,aes(y=value, x=x, color=variable, linetype=type)) +
  scale_y_continuous(name=unname(TeX("$L^{-1}(v)$")),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_x_continuous(name=TeX("$v$"),breaks=seq(-1,1,0.2),limits=c(-1,1)) +
  scale_color_manual(name=TeX("NB$(p_{i},r=3)$"),
                     labels=cols,values=c("black","blue","red","green")) +
  scale_linetype_manual(name="",
                        labels=lines,values=c("dashed","dotted")) +
  theme(legend.position=c(0.75,0.3),
        legend.title=element_text(size=20),
        legend.key.size = unit(3,"line"),
        legend.text=element_text(size=20,margin=margin(l=5,unit="pt")),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        panel.grid.major = element_line(linetype = 'solid',linewidth=0.25,colour = "gray")) +
  guides(color = guide_legend(order=1), linetype = guide_legend(order=0)) +
  geom_line(cex=0.5)
dev.off()




