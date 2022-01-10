
# Nice colours? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")


# Directed NONWEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
NonW_ST_matrix_rivers

# Spatiotemporal matrix 
NonW_ST_matrix_out_out
#NonW_ST_matrix_plots

# Spatiotemporal connectivity 
NonW_ST_connectivity_value
#NonW_ST_connectivity_plot

# Spatiotemporal Out.closenness 
NonW_ST_directed_Ocloseness_rivers
NonW_ST_directed_Oclo_mean <- list()
NonW_ST_directed_Oclo_sd <- list()
for (river in 1:length(NonW_ST_directed_Ocloseness_rivers)) {
  NonW_ST_directed_Oclo_mean[[river]] <- apply(NonW_ST_directed_Ocloseness_rivers[[river]],2,mean)
  NonW_ST_directed_Oclo_sd[[river]] <- apply(NonW_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

#NonW_ST_oclo_plot

# Spatiotemporal ALL.closenness 
NonW_ST_directed_Allcloseness_rivers
NonW_ST_directed_Allclo_mean <- list()
NonW_ST_directed_Allclo_sd <- list()
for (river in 1:length(NonW_ST_directed_Allcloseness_rivers)) {
  NonW_ST_directed_Allclo_mean[[river]] <- apply(NonW_ST_directed_Allcloseness_rivers[[river]],2,mean)
  NonW_ST_directed_Allclo_sd[[river]] <- apply(NonW_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

#NonW_ST_Allclo_plot

# Spatiotemporal Betweenness 
NonW_ST_directed_betweennes_rivers
NonW_ST_directed_betw_mean <- list()
NonW_ST_directed_betw_sd <- list()
for (river in 1:length(NonW_ST_directed_betweennes_rivers)) {
  NonW_ST_directed_betw_mean[[river]] <- apply(NonW_ST_directed_betweennes_rivers[[river]],2,mean)
  NonW_ST_directed_betw_sd[[river]] <- apply(NonW_ST_directed_betweennes_rivers[[river]],2,sd)
}

#NonW_ST_betw_plot


NonW_ST_directed_output <- data.frame(NonW_Dir_con=unlist(NonW_ST_connectivity_value),
                                      NonW_Dir_Ocl=unlist(NonW_ST_directed_Oclo_mean),
                                      NonW_Dir_Acl=unlist(NonW_ST_directed_Allclo_mean),
                                      NonW_Dir_Bet=unlist(NonW_ST_directed_betw_mean))


# Directed WEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
WEIG_ST_matrix_rivers 

# Spatiotemporal matrix 
WEIG_ST_matrix_out_out
#WEIG_ST_matrix_plots

# Spatiotemporal connectivity 
WEIG_ST_connectivity_value
#WEIG_ST_connectivity_plot

# Spatiotemporal Out.closenness 
WEIG_ST_directed_Ocloseness_rivers
WEIG_ST_Oclo_mean <- list()
WEIG_ST_Oclo_sd <- list()
for (river in 1:length(WEIG_ST_directed_Ocloseness_rivers)) {
  WEIG_ST_Oclo_mean[[river]] <- apply(WEIG_ST_directed_Ocloseness_rivers[[river]],2,mean)
  WEIG_ST_Oclo_sd[[river]] <- apply(WEIG_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

#WEIG_ST_Oclo_plot

# Spatiotemporal All.closenness 
WEIG_ST_directed_Allcloseness_rivers
WEIG_ST_Allclo_mean <- list()
WEIG_ST_Allclo_sd <- list()
for (river in 1:length(WEIG_ST_directed_Allcloseness_rivers)) {
  WEIG_ST_Allclo_mean[[river]] <- apply(WEIG_ST_directed_Allcloseness_rivers[[river]],2,mean)
  WEIG_ST_Allclo_sd[[river]] <- apply(WEIG_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

#WEIG_ST_Allclo_plot

# Spatiotemporal Betweenness 
WEIG_ST_directed_betweennes_rivers
WEIG_ST_betw_mean <- list()
WEIG_ST_betw_sd <- list()
for (river in 1:length(WEIG_ST_directed_betweennes_rivers)) {
  WEIG_ST_betw_mean[[river]] <- apply(WEIG_ST_directed_betweennes_rivers[[river]],2,mean)
  WEIG_ST_betw_sd[[river]] <- apply(WEIG_ST_directed_betweennes_rivers[[river]],2,sd)
}

#WEIG_ST_betw_plot


WEIG_ST_directed_output <- data.frame(WEIG_Dir_con=unlist(WEIG_ST_connectivity_value),
                                      WEIG_Dir_Ocl=unlist(WEIG_ST_Oclo_mean),
                                      WEIG_Dir_Acl=unlist(WEIG_ST_Allclo_mean),
                                      WEIG_Dir_Bet=unlist(WEIG_ST_betw_mean))


# Undirected river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
Un_NonW_ST_matrix_rivers

# Spatiotemporal matrix 
Un_NonW_ST_matrix_out_out
#Un_NonW_ST_matrix_plots

# Spatiotemporal matrix 
Un_NonW_ST_connectivity_value
#Un_NonW_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_NonW_ST_directed_Ocloseness_rivers
Un_NonW_ST_Oclo_mean <- list()
Un_NonW_ST_Oclo_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_Ocloseness_rivers)) {
  Un_NonW_ST_Oclo_mean[[river]] <- apply(Un_NonW_ST_directed_Ocloseness_rivers[[river]],2,mean)
  Un_NonW_ST_Oclo_sd[[river]] <- apply(Un_NonW_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

#Un_NonW_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_NonW_ST_directed_Allcloseness_rivers
Un_NonW_ST_Allclo_mean <- list()
Un_NonW_ST_Allclo_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_Allcloseness_rivers)) {
  Un_NonW_ST_Allclo_mean[[river]] <- apply(Un_NonW_ST_directed_Allcloseness_rivers[[river]],2,mean)
  Un_NonW_ST_Allclo_sd[[river]] <- apply(Un_NonW_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

#Un_NonW_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_NonW_ST_directed_betweennes_rivers
Un_NonW_ST_betw_mean <- list()
Un_NonW_ST_betw_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_betweennes_rivers)) {
  Un_NonW_ST_betw_mean[[river]] <- apply(Un_NonW_ST_directed_betweennes_rivers[[river]],2,mean)
  Un_NonW_ST_betw_sd[[river]] <- apply(Un_NonW_ST_directed_betweennes_rivers[[river]],2,sd)
}
#Un_NonW_ST_betw_plot

Un_NonW_ST_output <- data.frame(Un_NonW_con=unlist(Un_NonW_ST_connectivity_value),
                                Un_NonW_Ocl=unlist(Un_NonW_ST_Oclo_mean),
                                Un_NonW_Acl=unlist(Un_NonW_ST_Allclo_mean),
                                Un_NonW_Bet=unlist(Un_NonW_ST_betw_mean))


# Undirected WEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
Un_WEIG_ST_matrix_rivers

# Spatiotemporal matrix 
Un_WEIG_ST_matrix_out_out
#Un_WEIG_ST_matrix_plots

# Spatiotemporal matrix 
Un_WEIG_ST_connectivity_value
#Un_WEIG_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_WEIG_ST_directed_Ocloseness_rivers
Un_WEIG_ST_Oclo_mean <- list()
Un_WEIG_ST_Oclo_sd <- list()
for (river in 1:length(Un_WEIG_ST_directed_Ocloseness_rivers)) {
  Un_WEIG_ST_Oclo_mean[[river]] <- apply(Un_WEIG_ST_directed_Ocloseness_rivers[[river]],2,mean)
  Un_WEIG_ST_Oclo_sd[[river]] <- apply(Un_WEIG_ST_directed_Ocloseness_rivers[[river]],2,sd)
}
#Un_WEIG_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_WEIG_ST_directed_Allcloseness_rivers
Un_WEIG_ST_Allclo_mean <- list()
Un_WEIG_ST_Allclo_sd <- list()
for (river in 1:length(Un_WEIG_ST_directed_Allcloseness_rivers)) {
  Un_WEIG_ST_Allclo_mean[[river]] <- apply(Un_WEIG_ST_directed_Allcloseness_rivers[[river]],2,mean)
  Un_WEIG_ST_Allclo_sd[[river]] <- apply(Un_WEIG_ST_directed_Allcloseness_rivers[[river]],2,sd)
}
#Un_WEIG_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_WEIG_ST_directed_betweennes_rivers
Un_WEIG_ST_Betclo_mean <- list()
Un_WEIG_ST_Betclo_sd <- list()
for (river in 1:length(Un_WEIG_ST_directed_betweennes_rivers)) {
  Un_WEIG_ST_Betclo_mean[[river]] <- apply(Un_WEIG_ST_directed_betweennes_rivers[[river]],2,mean)
  Un_WEIG_ST_Betclo_sd[[river]] <- apply(Un_WEIG_ST_directed_betweennes_rivers[[river]],2,sd)
}
#Un_WEIG_ST_betw_plot 


Un_WEIG_ST_output <- data.frame(Un_WEIG_con=unlist(Un_WEIG_ST_connectivity_value),
                                Un_WEIG_Ocl=unlist(Un_WEIG_ST_Oclo_mean),
                                Un_WEIG_Acl=unlist(Un_WEIG_ST_Allclo_mean),
                                Un_WEIG_Bet=unlist(Un_WEIG_ST_Betclo_mean))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data processing and plotting __________________________________________________________________ ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Extract number of "HOBOS" per river and pre
n_each_riv <- c()
HOB_riv_ID <- c()
ups_dos <- c()
for (n in 1:length(Sites_list)) {
  # Vector with n sites 
  n_each_riv[n] <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv[n])
  ups_dos <- c(ups_dos,val_ups_dos)
  # ID of each stream
  ID <- Sites_list[[n]][,1]
  HOB_riv_ID <- c(HOB_riv_ID,ID)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Matrix HOBOs values__________________________________________________________________ ####
library(ape)
library(cowplot)

NonW_ST_directed_MatrixOut <- data.frame()
for (n in 1:length(NonW_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(NonW_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a)))[,1],
                          y=scores(metaMDS(as.dist(a)))[,2],
                          ID=Sites_list[[n]][,1],
                          UtoD=val_ups_dos)  
  NonW_ST_directed_MatrixOut <- rbind(NonW_ST_directed_MatrixOut,taula_out)
}


WEIG_ST_MatrixOut <- data.frame()
for (n in 1:length(WEIG_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(WEIG_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a)))[,1],
                          y=scores(metaMDS(as.dist(a)))[,2],
                          ID=Sites_list[[n]][,1],
                          UtoD=val_ups_dos)  
  WEIG_ST_MatrixOut <- rbind(WEIG_ST_MatrixOut,taula_out)
}



Un_NonW_ST_MatrixOut <- data.frame()
for (n in 1:length(Un_NonW_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(Un_NonW_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a)))[,1],
                          y=scores(metaMDS(as.dist(a)))[,2],
                          ID=Sites_list[[n]][,1],
                          UtoD=val_ups_dos)  
  Un_NonW_ST_MatrixOut <- rbind(Un_NonW_ST_MatrixOut,taula_out)
}



Un_WEIG_ST_MatrixOut <- data.frame()
for (n in 1:length(Un_WEIG_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(Un_WEIG_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a)))[,1],
                          y=scores(metaMDS(as.dist(a)))[,2],
                          ID=Sites_list[[n]][,1],
                          UtoD=val_ups_dos)  
  Un_WEIG_ST_MatrixOut <- rbind(Un_WEIG_ST_MatrixOut,taula_out)
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


legend_plots<- get_legend(ggplot(Un_WEIG_ST_MatrixOut)+
                            geom_point(aes(x=x, y=y, fill=ID,size=UtoD), shape=21,colour="grey10", alpha=0.5)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            scale_size(name="Upstream downstream position")+
                            theme_classic()+theme(legend.direction = "horizontal",
                                                  legend.box="horizontal"))

png(filename =paste("Figure/Data_treat/","Matrix_NMDS",".png"), 
    width = 750*6, height = 505*6, 
    units = "px",res = 300) 
grid.arrange(
ggplot(NonW_ST_directed_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=UtoD), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(WEIG_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=UtoD), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  scale_y_continuous(label=scientific_10)+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(Un_NonW_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=UtoD), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(Un_WEIG_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=UtoD), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  scale_y_continuous(label=scientific_10)+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

 legend_plots, nrow=5,ncol=1)
dev.off()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Individual HOBOs values__________________________________________________________________ ####

NonW_ST_directed_out <- data.frame(ID=HOB_riv_ID,UtoD=ups_dos,NonW_ST_directed_output)
WEIG_ST_directed_out <- data.frame(ID=HOB_riv_ID,UtoD=ups_dos,WEIG_ST_directed_output)
Un_NonW_ST_out <- data.frame(ID=HOB_riv_ID,UtoD=ups_dos,Un_NonW_ST_output)
Un_WEIG_ST_out <- data.frame(ID=HOB_riv_ID,UtoD=ups_dos,Un_WEIG_ST_output)

# Correlation for each approximation 
corrmorant::corrmorant(NonW_ST_directed_out,rescale = "by_range")
corrmorant::corrmorant(WEIG_ST_directed_out,rescale = "by_range")
corrmorant::corrmorant(Un_NonW_ST_out,rescale = "by_range")
corrmorant::corrmorant(Un_WEIG_ST_out,rescale = "by_range")

names_approxim <- c("Dir_NonW_ST", 
                    "Dir_WEIG_ST",
                    "Un_NonW_ST", 
                    "Un_WEIG_ST")
#PCA
library(ggfortify)
png(filename =paste("Figure/Data_treat/","PCAs ST ind. data",".png"), 
    width = 700*4, height = 600*4, 
    units = "px",res = 300) 
grid.arrange(
  autoplot(prcomp(NonW_ST_directed_out[,2:6], center = T, scale. = T),
         loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
                loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[1]),
  
  autoplot(prcomp(WEIG_ST_directed_out[,2:6], center = T, scale. = T),
           loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
           loadings.colour = 'red', loadings.label.colour="black")+
    geom_point(size=2, alpha=0.7,shape=21,
               color="grey30", aes(fill=WEIG_ST_directed_out[,1]))+
    scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
    theme_bw()+labs(title =names_approxim[2]), 
  
  autoplot(prcomp(Un_NonW_ST_out[,2:5], center = T, scale. = T),
           loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
           loadings.colour = 'red', loadings.label.colour="black")+
    geom_point(size=2, alpha=0.7,shape=21,
               color="grey30", aes(fill=Un_NonW_ST_out[,1]))+
    scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
    theme_bw()+labs(title =names_approxim[3]), 
  
  autoplot(prcomp(Un_WEIG_ST_out[,2:5], center = T, scale. = T),
           loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
           loadings.colour = 'red', loadings.label.colour="black")+
    geom_point(size=2, alpha=0.7,shape=21,
               color="grey30", aes(fill=Un_WEIG_ST_out[,1]))+
    scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
    theme_bw()+labs(title =names_approxim[4])
)
dev.off()
  

# Total correlation
corrmorant::corrmorant(cbind(NonW_ST_directed_out,
                             WEIG_ST_directed_out,
                             Un_NonW_ST_out,
                             Un_WEIG_ST_out),rescale = "by_range")
# ST connectivity
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,1],
                             WEIG_ST_directed_out[,1],
                             Un_NonW_ST_out[,1],
                             Un_WEIG_ST_out[,1]),rescale = "by_range")
# ST Oclo
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,2],
                             WEIG_ST_directed_out[,2],
                             Un_NonW_ST_out[,2],
                             Un_WEIG_ST_out[,2]),rescale = "by_range")

# ST Allclo
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,3],
                             WEIG_ST_directed_out[,3],
                             Un_NonW_ST_out[,3],
                             Un_WEIG_ST_out[,3]),rescale = "by_range")

# ST Bet
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,4],
                             WEIG_ST_directed_out[,4],
                             Un_NonW_ST_out[,4],
                             Un_WEIG_ST_out[,4]),rescale = "by_range")

four_approxim <- list(NonW_ST_directed_out, 
                      WEIG_ST_directed_out,
                      Un_NonW_ST_out, 
                      Un_WEIG_ST_out)

legend_plots<- get_legend(ggplot(NonW_ST_directed_out)+
                      geom_point(aes(x=UtoD,y=four_approxim[[a]][,3],fill=ID),shape=21, size=2)+
                      scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                      theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical"))

for (a in 1:length(four_approxim)) {
png(filename =paste("Figure/Data_treat/",names_approxim[a],".png"), 
      width = 700*4, height = 600*4, 
      units = "px",res = 300) 
grid.arrange(
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=UtoD,y=four_approxim[[a]][,3],fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=UtoD,y=four_approxim[[a]][,3],color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=UtoD,y=four_approxim[[a]][,3],color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Connecticity")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=UtoD,y=four_approxim[[a]][,4] ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=UtoD,y=four_approxim[[a]][,4] ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=UtoD,y=four_approxim[[a]][,4] ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Out closeness")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=UtoD,y=four_approxim[[a]][,5]  ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=UtoD,y=four_approxim[[a]][,5]  ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=UtoD,y=four_approxim[[a]][,5]  ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST All closeness")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=UtoD,y=four_approxim[[a]][,6]  ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=UtoD,y=four_approxim[[a]][,6]  ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=UtoD,y=four_approxim[[a]][,6]  ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Betwenness")+
    theme_classic()+theme(legend.position="none"),
  
  legend_plots, nrow=3 ,ncol=2 ,top=names_approxim[a])
dev.off()
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Final values__________________________________________________________________ ####
NonW_ST_matrix_out_out
WEIG_ST_matrix_out_out
Un_NonW_ST_matrix_out_out
Un_WEIG_ST_matrix_out_out

NonW_ST_directed_out
WEIG_ST_directed_out
Un_NonW_ST_out
Un_WEIG_ST_out



