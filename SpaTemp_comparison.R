
# Nice colors? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

## IMPORTANT MESSAGE: YOU NEED TO RUN the treatment script
# Needed values from the script "SpaTemp_HOBOS_treatment.R"
## RUN BEFORE PROCEEDING WITH THIS ONE

## You need to run the previous script in order to built the HOBOS dataset (a list with the HOBOS info for each river)
HOBOS_sites
# This two other objects are also built in the "SpaTemp_HOBOS_treatment.R" script. They contain the names 
## of each river and the order upstream to downstream within it. 
HOB_riv_ID
ups_dos

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Data processing and plotting for the different values 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Matrix ST outputs obtained from "SpaTemp_HOBOS_treatment.R"
NonW_ST_matrix_out_out
WEIG_ST_matrix_out_out
Un_NonW_ST_matrix_out_out
Un_WEIG_ST_matrix_out_out

# INDIVIDUAL ST outputs obtained from "SpaTemp_HOBOS_treatment.R"
NonW_ST_directed_out
WEIG_ST_directed_out
Un_NonW_ST_out
Un_WEIG_ST_out

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1. Matrix HOBOS values NMDS ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

library(ape)
library(vegan)
library(cowplot)

NonW_ST_directed_MatrixOut <- data.frame()
for (n in 1:length(NonW_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(NonW_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a),distance = "euclidian"))[,1],
                          y=scores(metaMDS(as.dist(a),distance = "euclidian"))[,2],
                          ID=Sites_list[[n]][,1],
                          DtoU=val_ups_dos)  
  NonW_ST_directed_MatrixOut <- rbind(NonW_ST_directed_MatrixOut,taula_out)
}


WEIG_ST_MatrixOut <- data.frame()
for (n in 1:length(WEIG_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(WEIG_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a),distance = "euclidian"))[,1],
                          y=scores(metaMDS(as.dist(a),distance = "euclidian"))[,2],
                          ID=Sites_list[[n]][,1],
                          DtoU=val_ups_dos)  
  WEIG_ST_MatrixOut <- rbind(WEIG_ST_MatrixOut,taula_out)
}



Un_NonW_ST_MatrixOut <- data.frame()
for (n in 1:length(Un_NonW_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(Un_NonW_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a),distance = "euclidian"))[,1],
                          y=scores(metaMDS(as.dist(a),distance = "euclidian"))[,2],
                          ID=Sites_list[[n]][,1],
                          DtoU=val_ups_dos)  
  Un_NonW_ST_MatrixOut <- rbind(Un_NonW_ST_MatrixOut,taula_out)
}



Un_WEIG_ST_MatrixOut <- data.frame()
for (n in 1:length(Un_WEIG_ST_matrix_out_out)) {
  n_each_riv <- nrow(Sites_list[[n]])
  # Sequence from upstream to downstream
  val_ups_dos <- seq(1,n_each_riv)
  
  a <- t(Un_WEIG_ST_matrix_out_out[[n]])
  diag(a) <- 0
  taula_out <- data.frame(x=scores(metaMDS(as.dist(a),distance = "euclidian"))[,1],
                          y=scores(metaMDS(as.dist(a),distance = "euclidian"))[,2],
                          ID=Sites_list[[n]][,1],
                          DtoU=val_ups_dos)  
  Un_WEIG_ST_MatrixOut <- rbind(Un_WEIG_ST_MatrixOut,taula_out)
}

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}


legend_plots<- get_legend(ggplot(Un_WEIG_ST_MatrixOut)+
                            geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            scale_size(name="Upstream downstream position")+
                            theme_classic()+theme(legend.direction = "horizontal",
                                                  legend.box="horizontal"))

png(filename =paste("Figure/Data_treat/","Matrix_NMDS",".png"), 
    width = 750*6, height = 505*6, 
    units = "px",res = 300) 
grid.arrange(
ggplot(NonW_ST_directed_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(WEIG_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  scale_y_continuous(label=scientific_10)+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(Un_NonW_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

ggplot(Un_WEIG_ST_MatrixOut)+
  geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
  stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI")+
  scale_color_CUNILLERA(palette = "LGTBI")+
  scale_y_continuous(label=scientific_10)+
  theme_classic()+theme(legend.position="none")+
  facet_grid(. ~ ID),

 legend_plots, nrow=5,ncol=1)
dev.off()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 2. Individual HOBOs values ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Correlation for each approximation 
png(filename ="Figure/Data_treat/Corrplots_Indiv.png", 
    width = 750*6, height = 700*6, 
    units = "px",res = 300) 
grid.arrange(
corrmorant::corrmorant(NonW_ST_directed_out,rescale = "by_range"),
corrmorant::corrmorant(WEIG_ST_directed_out,rescale = "by_range"),
corrmorant::corrmorant(Un_NonW_ST_out,rescale = "by_range"),
corrmorant::corrmorant(Un_WEIG_ST_out,rescale = "by_range"),
nrow=2,ncol=2)
dev.off()

names_approxim <- c("Dir_NonW_ST", 
                    "Dir_WEIG_ST",
                    "Un_NonW_ST", 
                    "Un_WEIG_ST")

#PCA Only between ST variables 
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
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,3:6],
                             WEIG_ST_directed_out[,3:6],
                             Un_NonW_ST_out[,3:6],
                             Un_WEIG_ST_out[,3:6]),rescale = "by_range")

png(filename ="Figure/Data_treat/Corrplots_Per_Index.png", 
    width = 750*6, height = 700*6, 
    units = "px",res = 300) 
grid.arrange(
# ST connectivity
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,3],
                             WEIG_ST_directed_out[,3],
                             Un_NonW_ST_out[,3],
                             Un_WEIG_ST_out[,3]),rescale = "by_range"),
# ST Oclo
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,4],
                             WEIG_ST_directed_out[,4],
                             Un_NonW_ST_out[,4],
                             Un_WEIG_ST_out[,4]),rescale = "by_range"),

# ST Allclo
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,5],
                             WEIG_ST_directed_out[,5],
                             Un_NonW_ST_out[,5],
                             Un_WEIG_ST_out[,5]),rescale = "by_range"),

# ST Bet
corrmorant::corrmorant(cbind(NonW_ST_directed_out[,6],
                             WEIG_ST_directed_out[,6],
                             Un_NonW_ST_out[,6],
                             Un_WEIG_ST_out[,6]),rescale = "by_range"),
nrow=2,ncol=2)
dev.off()


four_approxim <- list(NonW_ST_directed_out, 
                      WEIG_ST_directed_out,
                      Un_NonW_ST_out, 
                      Un_WEIG_ST_out)


for (a in 1:length(four_approxim)) {
png(filename =paste("Figure/Data_treat/",names_approxim[a],".png"), 
      width = 700*4, height = 600*4, 
      units = "px",res = 300)

rows_filter <- unlist(c(four_approxim[[a]]%>%
                          mutate(files=1:nrow(.))%>%
                          group_by(ID)%>%
                          mutate(maxim=max(DtoU), value=DtoU/maxim)%>%
                          filter(value==1)%>%
                          ungroup()%>%
                          dplyr::select(files)))
    
grid.arrange(
  ggplot(data = four_approxim[[a]][-rows_filter,])+
    geom_point(aes(x=DtoU,y=four_approxim[[a]][-rows_filter,3],fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=DtoU,y=four_approxim[[a]][-rows_filter,3],color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=DtoU,y=four_approxim[[a]][-rows_filter,3],color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Connecticity")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=DtoU,y=four_approxim[[a]][,4] ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=DtoU,y=four_approxim[[a]][,4] ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=DtoU,y=four_approxim[[a]][,4] ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Out closeness")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=DtoU,y=four_approxim[[a]][,5]  ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=DtoU,y=four_approxim[[a]][,5]  ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=DtoU,y=four_approxim[[a]][,5]  ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST All closeness")+
    theme_classic()+theme(legend.position="none"),
  
  ggplot(data = four_approxim[[a]])+
    geom_point(aes(x=DtoU,y=four_approxim[[a]][,6]  ,fill=ID),color="grey30", alpha=0.5, shape=21, size=2)+
    geom_line(aes(x=DtoU,y=four_approxim[[a]][,6]  ,color=ID),alpha=0.3, linetype=2)+
    geom_smooth(aes(x=DtoU,y=four_approxim[[a]][,6]  ,color=ID),alpha=0.2, method = "lm",se = F)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab("HOBO relative position Upstream-Downstream")+
    ylab("ST Betwenness")+
    theme_classic()+theme(legend.position="none"),
  
  get_legend(ggplot(NonW_ST_directed_out)+
               geom_point(aes(x=DtoU,y=four_approxim[[a]][,3],fill=ID),shape=21, size=2)+
               scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
               theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical")),
  nrow=3 ,ncol=2 ,top=names_approxim[a])
dev.off()
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 3. Comparing and Calculating OLD HOBOS values ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
source("Old_HOBOS_calculation.R")
Old_HOBOS_comb

# Printing image for a PCA with all the variables together 
for (a in 1:length(four_approxim)) {
  png(filename =paste("Figure/Data_treat/","PCAs ST ind. data vs TOT",".png"), 
      width = 700*4, height = 600*4, 
      units = "px",res = 300)
grid.arrange(
autoplot(prcomp(
         data.frame(NonW_ST_directed_out[,2:6],Old_HOBOS_comb[,3:5]),
         center = T, scale. = T),
         loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
         loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[1]),

autoplot(prcomp(
  data.frame(WEIG_ST_directed_out[,2:6],Old_HOBOS_comb[,3:5]),
  center = T, scale. = T),
  loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
  loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[2]),

autoplot(prcomp(
  data.frame(Un_NonW_ST_out[,2:5],Old_HOBOS_comb[,3:5]),
  center = T, scale. = T),
  loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
  loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[3]),

autoplot(prcomp(
  data.frame(Un_WEIG_ST_out[,2:5],Old_HOBOS_comb[,3:5]),
  center = T, scale. = T),
  loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
  loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[4])
)
dev.off()
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 4. Correlation plot with all values together to compare ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
png(filename =paste("Figure/Data_treat/CorPlot_STcon_OldHOB.png"), 
    width = 800*4, height = 600*4, 
    units = "px",res = 300)
corrmorant::corrmorant(cbind("STcon 1.1."=NonW_ST_directed_out$NonW_Dir_con,
                             "log(STcon) 1.2."=log(WEIG_ST_directed_out$WEIG_Dir_con+1),
                             "STcon 2.1"=Un_NonW_ST_out$Un_NonW_con,
                             "log(STcon) 2.2"=log(Un_WEIG_ST_out$Un_WEIG_con+1),
                             "log(TotDur)"= log(Old_HOBOS_comb[,3]+1),
                             "log(TotNum)"= log(Old_HOBOS_comb[,4]+1),
                             "log(TotLeng)"= log(Old_HOBOS_comb[,5]+1)),
                            rescale = "by_sd",corr_method = "spearman")
dev.off()

cor.test(NonW_ST_directed_out$NonW_Dir_con,
         log(Old_HOBOS_comb[,3]+1),method="spearman", exact = FALSE)
cor.test(NonW_ST_directed_out$NonW_Dir_con,
         log(Old_HOBOS_comb[,4]+1),method="spearman", exact = FALSE)

cor.test(Un_NonW_ST_out$Un_NonW_con,
         log(Old_HOBOS_comb[,3]+1),method="spearman",exact = FALSE)
cor.test(Un_NonW_ST_out$Un_NonW_con,
         log(Old_HOBOS_comb[,4]+1),method="spearman",exact = FALSE)
  




