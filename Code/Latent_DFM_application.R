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
Estim <- Latent_DFM_Estim(r=5,p,Cov_Z,identy_opt=2,shift=TRUE)


#################################################################################
# Forecasting:
#################################################################################
# Set up:
model <- list()
model$Lambda <- Estim$Lambda
model$Cov_eps <- diag(abs(Estim$Cov_eps))
model$Psi <- Estim$Psi
model$Cov_eta <- Estim$Cov_eta
model$Param <- Link_param$Param

#---------------------------------------------------#
# Forecasting:
Forecast <- Latent_DFM_Forecast(model,p,r=5,d,TT,X_t,H,dist,dist_opt)


#################################################################################
# Illustration:
#################################################################################
# Data:
reorder <- c(1,c(6:10),2,c(11:15),3,c(16:20),4,c(21:25),5,c(26:30))
label <- rownames(data)[reorder]

df_X <- data.frame(t(data[reorder,]),time=c(1:90))
df_melt <- melt(df_X,id="time")
df_melt$factor <- rep(c("C1","C2","C3","C4","C5"),each=90*6)
df_melt$factor <- factor(df_melt$factor,levels=c("C1","C2","C3","C4","C5"))
df_melt$variable <- factor(df_melt$variable,levels=label)
Mycol <- c(brewer.pal(n=7,name="Blues")[-1],
           brewer.pal(n=7,name="Greens")[-1],
           brewer.pal(n=7,name="Oranges")[-1],
           brewer.pal(n=7,name="Purples")[-1],
           brewer.pal(n=7,name="Greys")[-1])

# Plotting:
pdf(width=12,height=12,file="observation.pdf")
par(mar=c(0,0,0,0))
ggplot(df_melt,aes(x=time,y=value,color=variable)) +
  geom_line(aes(y=value,x=time)) +
  geom_point(aes(y=value,x=time)) +
  scale_x_continuous(name="Time",
                     breaks=c(1:90)[c(which(c(1:90)%%5==1),90)],
                     labels=c(1:90)[c(which(c(1:90)%%5==1),90)]) +
  scale_y_continuous(name=TeX("$\\X_{i,t}$"),limits=c(1,5)) + 
  scale_color_manual(name="Item",values=Mycol) +
  facet_grid(factor~.) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        axis.ticks.x=element_blank(),
        strip.text.x=element_text(size=9),strip.text.y=element_text(size=8))
dev.off()

#---------------------------------------------------#
# Estimation: Covariance matrices
df_cov_X0 <- data.frame(value=as.vector(Cov_X[reorder,reorder,(p+1)]),
                        row = rep(c(1:30),30),
                        column = rev(rep(c(1:30),1,each=30)))
df_cov_X1 <- data.frame(value=as.vector(Cov_X[reorder,reorder,(p+2)]),
                        row = rep(c(1:30),30),
                        column = rev(rep(c(1:30),1,each=30)))
df_cov_Z0 <- data.frame(value=as.vector(Cov_Z[reorder,reorder,(p+1)]),
                        row = rep(c(1:30),30),
                        column = rev(rep(c(1:30),1,each=30)))
df_cov_Z1 <- data.frame(value=as.vector(Cov_Z[reorder,reorder,(p+2)]),
                        row = rep(c(1:30),30),
                        column = rev(rep(c(1:30),1,each=30)))


pl_cov_X0 <- ggplot(df_cov_X0,aes(x=row,y=column,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Sigma}_{X}(0)$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="",
                   limits = rev(levels(as.factor(df_cov_X0$row))), 
                   breaks = levels(as.factor(df_cov_X0$row))[which(c(1:30)%%2==0)],
                   labels = label[which(c(1:30)%%2==0)])+
  scale_x_discrete(name="",
                   limits = levels(as.factor(df_cov_X0$column)), 
                   breaks = levels(as.factor(df_cov_X0$column))[which(c(1:30)%%2==1)],
                   labels = label[which(c(1:30)%%2==1)])+
  theme(axis.text.x=element_text(angle=30, vjust=0.7),
        panel.background = element_rect(fill = 'transparent'),
        legend.position="bottom")

pl_cov_X1 <- ggplot(df_cov_X1,aes(x=row,y=column,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Sigma}_{X}(1)$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="",
                   limits = rev(levels(as.factor(df_cov_X1$row))), 
                   breaks = levels(as.factor(df_cov_X1$row))[which(c(1:30)%%2==0)],
                   labels = label[which(c(1:30)%%2==0)])+
  scale_x_discrete(name="",
                   limits = levels(as.factor(df_cov_X1$column)), 
                   breaks = levels(as.factor(df_cov_X1$column))[which(c(1:30)%%2==1)],
                   labels = label[which(c(1:30)%%2==1)])+
  theme(axis.text.x=element_text(angle=30, vjust=0.7),
        panel.background = element_rect(fill = 'transparent'),
        legend.position="bottom")

pl_cov_Z0 <- ggplot(df_cov_Z0,aes(x=row,y=column,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Sigma}_{Z}(0)
              =\\widehat{L}^{-1}(\\widehat{\\Sigma}_{X}(0))$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="",
                   limits = rev(levels(as.factor(df_cov_Z0$row))), 
                   breaks = levels(as.factor(df_cov_Z0$row))[which(c(1:30)%%2==0)],
                   labels = label[which(c(1:30)%%2==0)])+
  scale_x_discrete(name="",
                   limits = levels(as.factor(df_cov_Z0$column)), 
                   breaks = levels(as.factor(df_cov_Z0$column))[which(c(1:30)%%2==1)],
                   labels = label[which(c(1:30)%%2==1)])+
  theme(axis.text.x=element_text(angle=30, vjust=0.7),
        panel.background = element_rect(fill = 'transparent'),
        legend.position="bottom")

pl_cov_Z1 <- ggplot(df_cov_Z1,aes(x=row,y=column,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Sigma}_{Z}(1)
              =\\widehat{L}^{-1}(\\widehat{\\Sigma}_{X}(1))$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="",
                   limits = rev(levels(as.factor(df_cov_Z1$row))), 
                   breaks = levels(as.factor(df_cov_Z1$row))[which(c(1:30)%%2==0)],
                   labels = label[which(c(1:30)%%2==0)])+
  scale_x_discrete(name="",
                   limits = levels(as.factor(df_cov_Z1$column)), 
                   breaks = levels(as.factor(df_cov_Z1$column))[which(c(1:30)%%2==1)],
                   labels = label[which(c(1:30)%%2==1)])+
  theme(axis.text.x=element_text(angle=30, vjust=0.7),
        panel.background = element_rect(fill = 'transparent'),
        legend.position="bottom")

# Plotting:
pdf(width=12,height=12,file="covariances.pdf")
par(mar=c(0,0,0,0))
ggarrange(pl_cov_X0,pl_cov_X1,pl_cov_Z0,pl_cov_Z1,
          ncol=2,nrow=2,legend="bottom",common.legend=TRUE)
dev.off()

#---------------------------------------------------#
# Estimation: Loadings matrix and transition matrix
df_lambda <- data.frame(value = as.vector(model$Lambda[reorder,]),
                        item = rev(rep(c(1:30),5)),
                        factor = rep(c(1:5),5,each=30))
df_lambda$item <- factor(df_lambda$item,levels=c(1:30))
df_lambda$factor <- factor(df_lambda$factor,levels=c(1:5))

pl_lambda <- ggplot(df_lambda,aes(x=factor,y=item,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Lambda}$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="Item",
                   limits = levels(as.factor(df_lambda$item)), 
                   breaks = levels(as.factor(df_lambda$item)),
                   labels = label)+
  scale_x_discrete(name="Factor",
                   limits = levels(as.factor(df_lambda$factor)), 
                   breaks = levels(as.factor(df_lambda$factor)),
                   labels = c("C1","C2","C3","C4","C5")) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'))

df_psi <- data.frame(value = as.vector(model$Psi),
                     to = rep(c(1:5),5),
                     from = rep(c(1:5),1,each=5))
df_psi$from <- factor(df_psi$from,levels=c(1,2,3,4,5))
df_psi$to <- factor(df_psi$to,levels=c(1,2,3,4,5))

pl_psi <- ggplot(df_psi,aes(x=from,y=to,fill=value))+
  geom_tile()+
  ggtitle(TeX("$\\widehat{\\Psi}_{1}$"))+
  scale_fill_gradient2(name="",low="darkblue",high="darkred")+
  scale_y_discrete(name="",
                   limits = rev(levels(as.factor(df_psi$to))), 
                   breaks= levels(as.factor(df_psi$to)),
                   labels = c("C1","C2","C3","C4","C5"))+
  scale_x_discrete(name="",
                   limits = levels(as.factor(df_psi$from)), 
                   breaks = levels(as.factor(df_psi$from)),
                   labels = c("C1","C2","C3","C4","C5")) +
  theme(panel.background = element_rect(fill = 'transparent'))


# Plotting:
pdf(width=12,height=8,file="parameters.pdf")
par(mar=c(0,0,0,0))
ggarrange(pl_lambda,pl_psi,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)
dev.off()

#---------------------------------------------------# 
# Forecasting: Particles
list_Z_T <- Forecast$Z_t[reorder,,]
list_Z_hat <- Forecast$Z_hat_h[reorder,,]

particle_T <- data.frame(value=as.vector(list_Z_T),
                         item=rep(rep(label,5),100),
                         time=rep(c(81:85),100,each=30),
                         factor=rep(rep(rep(c("C1","C2","C3","C4","C5"),each=6),10),100),
                         particle=rep(c(1:100),1,each=300))
particle_T$item <- factor(particle_T$item,levels=unique(particle_T$item))
particle_T$factor <- factor(particle_T$factor,levels=unique(particle_T$factor))

particle_T_hat <- data.frame(value=as.vector(list_Z_hat),
                         item=rep(rep(label,5),100),
                         time=rep(c(86:90),100,each=30),
                         factor=rep(rep(rep(c("C1","C2","C3","C4","C5"),each=6),10),100),
                         particle=rep(c(1:100),1,each=300))
particle_T_hat$item <- factor(particle_T_hat$item,levels=unique(particle_T_hat$item))
particle_T_hat$factor <- factor(particle_T_hat$factor,levels=unique(particle_T_hat$factor))

aug_thrs <- melt(t(mapply(x=1:30,function(x)Forecast$threshold[[x]]))[reorder,])[,-1]
aug_thrs$item <- rep(label,each=1,4)
aug_thrs$item <- factor(aug_thrs$item,levels=label)
aug_thrs$factor <- rep(rep(c("C1","C2","C3","C4","C5"),each=6),4)
aug_thrs$Var2 <- factor(aug_thrs$Var2,levels=unique(aug_thrs$Var2))

bdd <- min(max(abs(floor(min(list_Z_T))),
           abs(ceiling(max(list_Z_T))),
           abs(floor(min(list_Z_hat))),
           abs(ceiling(max(list_Z_hat)))),6)


pl_c1_T <- ggplot(particle_T[which(particle_T$factor=="C1"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Last 5 observation period ($t$)"),breaks=c(81:85)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,t}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[1:6]) +
  facet_grid(factor(item,levels=label[1:6])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C1"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))
  
pl_c1_T_hat <- ggplot(particle_T_hat[which(particle_T$factor=="C1"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Next 5 forecast period ($t$)"),breaks=c(86:90)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,T+h|T}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[1:6]) +
  facet_grid(factor(item,levels=label[1:6])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C1"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))
  
ggarrange(pl_c1_T,pl_c1_T_hat,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)

pl_c2_T <- ggplot(particle_T[which(particle_T$factor=="C2"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Last 5 observation period ($t$)"),breaks=c(81:85)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,t}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[7:12]) +
  facet_grid(factor(item,levels=label[7:12])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C2"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

pl_c2_T_hat <- ggplot(particle_T_hat[which(particle_T$factor=="C2"),],
                      aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Next 5 forecast period ($t$)"),breaks=c(86:90)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,T+h|T}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[7:12]) +
  facet_grid(factor(item,levels=label[7:12])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C2"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

ggarrange(pl_c2_T,pl_c2_T_hat,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)


pl_c3_T <- ggplot(particle_T[which(particle_T$factor=="C3"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Last 5 observation period ($t$)"),breaks=c(81:85)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,t}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[13:18]) +
  facet_grid(factor(item,levels=label[13:18])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C3"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

pl_c3_T_hat <- ggplot(particle_T_hat[which(particle_T$factor=="C3"),],
                      aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Next 5 forecast period ($t$)"),breaks=c(86:90)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,T+h|T}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[13:18]) +
  facet_grid(factor(item,levels=label[13:18])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C3"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

ggarrange(pl_c3_T,pl_c3_T_hat,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)


pl_c4_T <- ggplot(particle_T[which(particle_T$factor=="C4"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Last 5 observation period ($t$)"),breaks=c(81:85)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,t}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[19:24]) +
  facet_grid(factor(item,levels=label[19:24])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C4"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

pl_c4_T_hat <- ggplot(particle_T_hat[which(particle_T$factor=="C4"),],
                      aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Next 5 forecast period ($t$)"),breaks=c(86:90)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,T+h|T}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[19:24]) +
  facet_grid(factor(item,levels=label[19:24])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C4"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

ggarrange(pl_c4_T,pl_c4_T_hat,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)



pl_c5_T <- ggplot(particle_T[which(particle_T$factor=="C5"),],
                  aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Last 5 observation period ($t$)"),breaks=c(81:85)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,t}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[25:30]) +
  facet_grid(factor(item,levels=label[25:30])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C5"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

pl_c5_T_hat <- ggplot(particle_T_hat[which(particle_T$factor=="C5"),],
                      aes(x=time,y=value,color=item)) +
  geom_line(aes(group=particle),alpha=0.5,linewidth=0.1) +
  scale_x_continuous(name=TeX("Next 5 forecast period ($t$)"),breaks=c(86:90)) +
  scale_y_continuous(name=TeX("$\\hat{Z}_{i,T+h|T}^{(k)}$"),limits=c(-bdd,bdd)) +
  scale_color_manual(name="Item",values=Mycol[25:30]) +
  facet_grid(factor(item,levels=label[25:30])~factor) +
  geom_hline(aes(yintercept=value,group=item),linetype="dotted",
             data=aug_thrs[which(aug_thrs$factor=="C5"),],linewidth=0.5)+
  theme(legend.position="bottom",
        panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x=element_blank(),strip.text.y=element_blank()) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1),nrow=1))

ggarrange(pl_c5_T,pl_c5_T_hat,
          ncol=2,nrow=1,legend="bottom",common.legend=TRUE)


# Plotting:
pdf(width=13,height=13,file="particles.pdf")
par(mar=c(0,0,0,0))
ggarrange(ggarrange(pl_c1_T,pl_c1_T_hat,ncol=2,nrow=1,legend="bottom",common.legend=TRUE),
          ggarrange(pl_c2_T,pl_c2_T_hat,ncol=2,nrow=1,legend="bottom",common.legend=TRUE),
          ggarrange(pl_c3_T,pl_c3_T_hat,ncol=2,nrow=1,legend="bottom",common.legend=TRUE),
          ggarrange(pl_c4_T,pl_c4_T_hat,ncol=2,nrow=1,legend="bottom",common.legend=TRUE),
          ggarrange(pl_c5_T,pl_c5_T_hat,ncol=2,nrow=1,legend="bottom",common.legend=TRUE),
          nrow=3,ncol=2)
dev.off()


#---------------------------------------------------# 
# Forecasting: Evaluation
H <- 5
TT <- 85

X_ref <- data[,(TT+1):(TT+H)]

diff_sis <- abs(Forecast$X_t_h - X_ref)
diff_L <- abs(data[,TT] - X_ref)
M <- as.numeric(as.character(mapply(x=1:d,function(x)data.frame(sort(table(X_t[x,]),
                                                                     decreasing=TRUE))[1,1])))
diff_M <- abs(M - X_ref)

df_sis <- data.frame(value=as.vector(diff_sis[reorder,]),
                     item=rep(label,5),
                     horizon=rep(1:5,each=30),
                     method="SIS/R")
df_L <- data.frame(value=as.vector(diff_L[reorder,]),
                    item=rep(label,5),
                    horizon=rep(1:5,each=30),
                    method="Last")
df_M <- data.frame(value=as.vector(diff_M[reorder,]),
                    item=rep(label,5),
                    horizon=rep(1:5,each=30),
                    method="Marginal")
df_forecast <- rbind(df_sis,df_L,df_M)
df_forecast$method <- factor(df_forecast$method,levels=unique(df_forecast$method))
df_forecast$item <- factor(df_forecast$item,levels=unique(df_forecast$item))



# Plotting:
pdf(width=12,height=8,file="forecasting.pdf")
par(mar=c(0,0,0,0))
ggplot(df_forecast, aes(y=value,x=horizon,color=item),fill=method) + 
  theme(legend.position="bottom",
        axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        strip.text.x=element_text(size=9),strip.text.y=element_text(size=8)) +
  scale_y_continuous("Absolute difference",breaks=c(0:4),labels=c(0:4)) +
  scale_color_manual(name="Item",values=Mycol) +
  scale_x_discrete(name="5-step-ahead horizon") +
  geom_line(aes(y=value,x=horizon)) +
  geom_point(aes(y=value,x=horizon)) +
  facet_grid(method~item) + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.border = element_rect(fill=NA,color="black",linetype="solid"),
        strip.text.x = element_blank())
dev.off()



