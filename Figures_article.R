detach("package:plyr", unload = TRUE)
# Figure 2  ####

Sites_plot <- data.frame(Riera=unique(Sites$Riera)) %>% arrange(Riera) %>% filter(Riera!="A") %>% filter(Riera!="G")
for (site_pl in 1:nrow(Sites_plot)) {
png(filename =paste("Figure/Data_treat/","Figure","_",c(Sites_plot[,1])[site_pl],".png",sep =""), 
    width = 530*5.5, height = 670*5.5, 
    units = "px",res = 300)

Ste_cha <- theme(axis.text = element_blank(),
                         axis.ticks = element_blank(),
                         axis.line = element_blank(),
                         panel.border = element_rect(size=1, colour = "grey20", fill=NA))
HOBO3_highlight <- NULL

if(site_pl==7){
HOBO3_highlight <- geom_nodes(size=c(0,0,10,rep(0,7)),colour="red",alpha=c(0,0,0.5,rep(0,7)))
}

grid.arrange(
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$STcon_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="A) STcon",subtitle="DirBin")),
  arrangeGrob(Dir_WEIG_Net$STcon_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="DirWei")),
  arrangeGrob(UnD_NonW_Net$STcon_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="UndBin")),
  arrangeGrob(UnD_WEIG_Net$STcon_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="UndWei")),
  top=""
),
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$STconmat_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="B) STconmat",subtitle="DirBin")),
  arrangeGrob(Dir_WEIG_Net$STconmat_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="DirWei")),
  arrangeGrob(UnD_NonW_Net$STconmat_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="UndBin")),
  arrangeGrob(UnD_WEIG_Net$STconmat_plo[[site_pl]]+Ste_cha+HOBO3_highlight+ labs(title="",subtitle="UndWei")),
  top=""
))
dev.off()
}

#

library(plotly)
library(ggfortify)

# Figure 3 ####
data_PCA <- data.frame(data.frame("DirBin"=NonW_ST_directed_out$NonW_Dir_con,
                                  "DirWei"=WEIG_ST_directed_out$WEIG_Dir_con,
                                  "UndBin"=Un_NonW_ST_out$Un_NonW_con, 
                                  "UndWei"=Un_WEIG_ST_out$Un_WEIG_con),
            Old_HOBOS_comb[,3:5])


png(filename =paste("Figure/Data_treat/","Figure3",".png"), 
    width = 700*2.5, height = 550*2.5, 
    units = "px",res = 300)

PCA_info <- prcomp(data_PCA,center = T, scale. = T)
a <- data.frame(
     adonis2(PCA_info$x[,2]~NonW_ST_directed_out$ID*NonW_ST_directed_out$DtoU,
        method = "euclidean")
      )
rownames(a) <- c("Stream ID","Upstr-Downstr", "Interaction","Residual","Total")
caption_text <- paste(rownames(a)[1],"pvalue=",a[1,5]," ",
                       rownames(a)[2],"pvalue=",a[2,5]," ",
                       rownames(a)[3],"pvalue=",a[3,5], sep = " ")
autoplot(prcomp(data_PCA,
  center = T, scale. = T),
  loadings=T, loadings.label = TRUE, loadings.label.size = 6, shape = F, label=F,
  loadings.colour = 'grey15', loadings.label.colour="black")+
  geom_point(alpha=0.5,shape=21,
             color="grey30", 
             aes(size=NonW_ST_directed_out$DtoU ,fill=NonW_ST_directed_out$ID))+
  stat_ellipse(aes(colour=NonW_ST_directed_out$ID),type = "t", linetype=2,size=1, alpha=0.5)+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  scale_color_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  scale_size(name="Upstr-Downstr")+
  scale_x_continuous(limits = c(-0.35,0.4))+
  scale_y_continuous(limits = c(-0.32,0.32))+
  labs(caption = caption_text)+
  theme_bw()
dev.off()
#

# Figure 4  ####
png(filename =paste("Figure/Data_treat/Figure4.png"), 
      width = 700*6, height = 650*6, 
      units = "px",res = 350)
  rows_filter <- unlist(c(four_approxim[[1]]%>%
                            mutate(files=1:nrow(.))%>%
                            group_by(ID)%>%
                            mutate(maxim=max(DtoU), value=DtoU/maxim)%>%
                            filter(value==1)%>%
                            ungroup()%>%
                            dplyr::select(files)))
  
  alpha_stream <- as.numeric(unlist(c(four_approxim[[1]]%>%
                  mutate(alpha_stream=ifelse(ID=="VH", 0.8, 0.8))%>%
                  ungroup()%>%
                  dplyr::select(alpha_stream))))
  alpha_HOBOS3 <- c(rep(0,58),1,rep(0,7))
  
  grid.arrange(
  #A1 line plot  
  arrangeGrob(
    ggplot(data = four_approxim[[1]][-rows_filter,])+
      #geom_point(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3]),color="red",alpha=alpha_HOBOS3[-rows_filter], shape=16, size=5,)+
      geom_point(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],fill=ID),color="grey30", alpha=alpha_stream[-rows_filter], shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],color=ID),alpha=alpha_stream[-rows_filter], linetype=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],color=ID, alpha=ID), 
                size=1.3,stat="smooth",method = "lm")+
      scale_alpha_manual(values = c(rep(0.3,6),0.9))+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "DirBin - STcon", title = "A)")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      theme_classic()+theme(legend.position="none"),
    top = ""),
  
  #A1 NMDS plot  
  arrangeGrob(    
    ggplot(NonW_ST_directed_MatrixOut)+
      #geom_point(aes(x=x,y=y, fill=ID,size=DtoU),color="red",alpha=alpha_HOBOS3, shape=16, size=5,)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "DirBin - STconmat", title = "B)")+
      scale_y_continuous(label=NULL)+
      scale_x_continuous(label=NULL)+
      theme_classic()+theme(legend.position="none",axis.ticks = element_blank())+
      facet_grid(. ~ ID),
    top = ""),
    
  
  #A2 line plot  
  arrangeGrob( 
    ggplot(data = four_approxim[[2]][,])+
      #geom_point(aes(x=DtoU,y=four_approxim[[2]][,3]),color="red",alpha=alpha_HOBOS3, shape=16, size=5,)+
      geom_point(aes(x=DtoU,y=four_approxim[[2]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[2]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[2]][,3],color=ID, alpha=ID), 
                size=1.3,stat="smooth",method = "lm")+
      scale_alpha_manual(values = c(rep(0.3,6),0.9))+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "DirWei - STcon")+ 
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      scale_y_continuous(limits = c(0,15000))+
      theme_classic()+theme(legend.position="none"),
    top=""),
    
  #A2 NMDS plot
  arrangeGrob(
    ggplot(WEIG_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=NULL)+
      scale_x_continuous(label=NULL)+
      labs(subtitle = "DirWei - STconmat")+
      theme_classic()+theme(legend.position="none",axis.ticks = element_blank())+
    facet_grid(. ~ ID), top = ""),
  
  #B1 line plot
  arrangeGrob(
    ggplot(data = four_approxim[[3]][,])+
      #geom_point(aes(x=DtoU,y=four_approxim[[3]][,3]),color="red",alpha=alpha_HOBOS3, shape=16, size=5,)+
      geom_point(aes(x=DtoU,y=four_approxim[[3]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[3]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[3]][,3],color=ID, alpha=ID), 
                size=1.3,stat="smooth",method = "lm")+
      scale_alpha_manual(values = c(rep(0.3,6),0.9))+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "UndBin - STcon")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      scale_y_continuous(limits = c(0,11.5))+
      theme_classic()+theme(legend.position="none"),
    top=""),
  
  #B1 NMDS plot
  arrangeGrob(
    ggplot(Un_NonW_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "UndBin - STconmat")+
      scale_y_continuous(label=NULL)+
      scale_x_continuous(label=NULL)+
      theme_classic()+theme(legend.position="none",axis.ticks = element_blank())+
      facet_grid(. ~ ID),top = ""),
  
  #B2 line plot  
  arrangeGrob(
    ggplot(data = four_approxim[[4]][,])+
      #geom_point(aes(x=DtoU,y=four_approxim[[4]][,3]),color="red",alpha=alpha_HOBOS3, shape=16, size=5,)+
      geom_point(aes(x=DtoU,y=four_approxim[[4]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[4]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[4]][,3],color=ID, alpha=ID), 
                size=1.3,stat="smooth",method = "lm")+
      scale_alpha_manual(values = c(rep(0.3,6),0.9))+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      labs(subtitle = "UndWei - STcon")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      scale_y_continuous(limits = c(0,60000))+
      theme_classic()+theme(legend.position="none"),
    top=""),
  
  #B1 NMDS plot
  arrangeGrob(
    ggplot(Un_WEIG_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=NULL)+
      scale_x_continuous(label=NULL)+
      labs(subtitle = "UndWei - STconmat")+
      theme_classic()+theme(legend.position="none",axis.ticks = element_blank())+
      facet_grid(. ~ ID),top = ""),
  

    get_legend(ggplot(NonW_ST_directed_out)+
                 geom_point(aes(x=DtoU,y=four_approxim[[1]][,3],fill=ID),shape=21, size=5)+
                 scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                 theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical")),
    nrow=5 ,ncol=2,widths=c(0.6,1.1))
  dev.off()

#Figure 5  ####  

# Select each significant plot.   
plots_BID_sign <- list(# Plots sign richness
                       plots_HOB_BDD_total[[1]][[3]]+labs(title="", subtitle = "Scenario UndBin Active dispersers" ),
                       plots_HOB_BDD_total[[1]][[6]]+labs(title="A) Richness", subtitle = "Scenario DirWei Passive dispersers"),
                       plots_HOB_BDD_total[[1]][[7]]+labs(title="", subtitle = "Scenario UndBin Passive dispersers"),
                       # Shannon
                       plots_HOB_BDD_total[[2]][[3]]+labs(title="", subtitle = "Scenario UndBin Active dispersers" ),
                       plots_HOB_BDD_total[[2]][[7]]+labs(title="B) Shannon", subtitle = "Scenario UndBin Passive dispersers" ),
                       # Trait ab.
                       plots_HOB_BDD_total[[3]][[7]]+labs(title="C) Trait abundance", subtitle = "Scenario UndBin Passive dispersers"),
                       #BrayC
                       plots_HOB_BDD_total[[4]][[2]]+labs(title="", subtitle = "Scenario DirWei Active dispersers"),
                       plots_HOB_BDD_total[[4]][[4]]+labs(title="", subtitle = "Scenario UndWei Active dispersers"),
                       plots_HOB_BDD_total[[4]][[6]]+labs(title="D) Bray-curtis", subtitle = "Scenario DirWei Passive dispersers"),
                       plots_HOB_BDD_total[[4]][[8]]+labs(title="", subtitle = "Scenario UndWei Passive dispersers"),
                       #Jaccard
                       plots_HOB_BDD_total[[5]][[1]]+labs(title="", subtitle = "Scenario DirBin Active dispersers"),
                       plots_HOB_BDD_total[[5]][[4]]+labs(title="", subtitle = "Scenario UndWei Active dispersers"),
                       plots_HOB_BDD_total[[5]][[5]]+labs(title="E) Jaccard", subtitle = "Scenario DirBin Passive dispersers"),
                       plots_HOB_BDD_total[[5]][[7]]+labs(title="", subtitle = "Scenario UndBin Passive dispersers"))
                      
legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "horizontal",
                                                  legend.box="horizontal"))


png(filename =paste("Figure/Data_treat/Figure5.png"), 
    width = 850*4.5, height = 750*4.5, 
    units = "px",res = 300)
grid.arrange(
  arrangeGrob( 
    plots_BID_sign[[2]], 
    plots_BID_sign[[3]],
    plots_BID_sign[[1]],
    legend_plots,
    ncol=4, top=""), 
  
  arrangeGrob( 
    plots_BID_sign[[5]],
    plots_BID_sign[[4]],
    plots_BID_sign[[6]],
    ncol=4, top=""),
  
  arrangeGrob( 
    plots_BID_sign[[9]],
    plots_BID_sign[[7]],
    plots_BID_sign[[10]],
    plots_BID_sign[[8]],
    ncol=4, top=""),
  
  arrangeGrob( 
    plots_BID_sign[[13]],
    plots_BID_sign[[11]],
    plots_BID_sign[[14]],
    plots_BID_sign[[12]],
    ncol=4, top=""),

  nrow=4)
dev.off()

