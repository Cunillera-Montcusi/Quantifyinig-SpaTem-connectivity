

png(filename =paste("Figure/Data_treat/","Fig2",".png"), 
    width = 530*5, height = 670*5, 
    units = "px",res = 300)
grid.arrange(
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$Dir_NonW_ST_con_plo[[7]],bottom="A1) Directed Binary"),
  arrangeGrob(Dir_WEIG_Net$Dir_WEIG_ST_con_plo[[7]],bottom="A2) Directed Weighted"),
  arrangeGrob(UnD_NonW_Net$UnD_NonW_ST_con_plo[[7]],bottom="B1) Undirected Binary"),
  arrangeGrob(UnD_WEIG_Net$UnD_WEIG_ST_con_plo[[7]],bottom="B2) Undirected Weighted"),
  top="A) STcon"
),
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$Dir_NonW_ST_matrix_plo[[7]],bottom="A1) Directed Binary"),
  arrangeGrob(Dir_WEIG_Net$Dir_WEIG_ST_matrix_plo[[7]],bottom="A2) Directed Weighted"),
  arrangeGrob(UnD_NonW_Net$UnD_NonW_ST_matrix_plo[[7]],bottom="B1) Undirected Binary"),
  arrangeGrob(UnD_WEIG_Net$UnD_WEIG_ST_matrix_plo[[7]],bottom="B2) Undirected Weighted"),
  top="B) STconmat"
))
dev.off()


data_PCA <- data.frame(data.frame("STcon A1."=NonW_ST_directed_out$NonW_Dir_con,
                                  "STcon A2."=WEIG_ST_directed_out$WEIG_Dir_con,
                                  "STcon B1."=Un_NonW_ST_out$Un_NonW_con, 
                                  "STcon B2."=Un_WEIG_ST_out$Un_WEIG_con),
            Old_HOBOS_comb[,3:5])


png(filename =paste("Figure/Data_treat/","PCAs_Fig3",".png"), 
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



png(filename =paste("Figure/Data_treat/Figure4.png"), 
      width = 700*4, height = 600*4, 
      units = "px",res = 300)
  rows_filter <- unlist(c(four_approxim[[1]]%>%
                            mutate(files=1:nrow(.))%>%
                            group_by(ID)%>%
                            mutate(maxim=max(DtoU), value=DtoU/maxim)%>%
                            filter(value==1)%>%
                            ungroup()%>%
                            select(files)))
  
  alpha_stream <- as.numeric(unlist(c(four_approxim[[1]]%>%
                  mutate(alpha_stream=ifelse(ID=="VH", 1, 0.3))%>%
                  ungroup()%>%
                  select(alpha_stream))))

  
  grid.arrange(
  arrangeGrob(
    ggplot(data = four_approxim[[1]][-rows_filter,])+
      geom_point(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],fill=ID),color="grey30", alpha=alpha_stream[-rows_filter], shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],color=ID),alpha=alpha_stream[-rows_filter], linetype=2)+
      geom_smooth(aes(x=DtoU,y=four_approxim[[1]][-rows_filter,3],color=ID),alpha=0.6, method = "lm",se = F)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      theme_classic()+theme(legend.position="none"),
    top="1.1) Directed Binary"),
    
  
  arrangeGrob(
    ggplot(data = four_approxim[[2]][,])+
      geom_point(aes(x=DtoU,y=four_approxim[[2]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[2]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_smooth(aes(x=DtoU,y=four_approxim[[2]][,3],color=ID),alpha=0.6, method = "lm",se = F)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      theme_classic()+theme(legend.position="none"),
  top="1.2) Directed Weighted"),
  
  arrangeGrob(
    ggplot(data = four_approxim[[3]][,])+
      geom_point(aes(x=DtoU,y=four_approxim[[3]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[3]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_smooth(aes(x=DtoU,y=four_approxim[[3]][,3],color=ID),alpha=0.6, method = "lm",se = F)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      theme_classic()+theme(legend.position="none"),
  top="2.1) Undirected Binary"),
    
  arrangeGrob(
    ggplot(data = four_approxim[[4]][,])+
      geom_point(aes(x=DtoU,y=four_approxim[[4]][,3],fill=ID),color="grey30", alpha=alpha_stream, shape=21, size=2)+
      geom_line(aes(x=DtoU,y=four_approxim[[4]][,3],color=ID),alpha=alpha_stream, linetype=2)+
      geom_smooth(aes(x=DtoU,y=four_approxim[[4]][,3],color=ID),alpha=0.6, method = "lm",se = F)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      xlab("HOBO relative position Upstream-Downstream")+
      ylab("ST Connecticity")+
      theme_classic()+theme(legend.position="none"),
    top="2.2) Undirected Weighted"),
    
    get_legend(ggplot(NonW_ST_directed_out)+
                 geom_point(aes(x=DtoU,y=four_approxim[[1]][,3],fill=ID),shape=21, size=2)+
                 scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                 theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical")),
    nrow=3 ,ncol=2)
  dev.off()

  legend_plots<- get_legend(ggplot(Un_WEIG_ST_MatrixOut)+
                              geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=0.5)+
                              scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                              scale_size(name="Upstream downstream position")+
                              theme_classic()+theme(legend.direction = "horizontal",
                                                    legend.box="vertical"))
  
  png(filename =paste("Figure/Data_treat/Figure4.png"), 
      width = 650*6, height = 505*6, 
      units = "px",res = 300) 
  grid.arrange(
    ggplot(NonW_ST_directed_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=scientific_10)+
      scale_x_continuous(label=scientific_10)+
      theme_classic()+theme(legend.position="none"),
    
    ggplot(WEIG_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=scientific_10)+
      scale_x_continuous(label=scientific_10)+
      theme_classic()+theme(legend.position="none"),
    
    ggplot(Un_NonW_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=scientific_10)+
      scale_x_continuous(label=scientific_10)+
      theme_classic()+theme(legend.position="none"),
    
    ggplot(Un_WEIG_ST_MatrixOut)+
      geom_point(aes(x=x, y=y, fill=ID,size=DtoU), shape=21,colour="grey10", alpha=alpha_stream)+
      stat_ellipse(aes(x=x, y=y, colour=ID),type = "t", linetype=2,size=1, alpha=0.5)+
      scale_fill_CUNILLERA(palette = "LGTBI")+
      scale_color_CUNILLERA(palette = "LGTBI")+
      scale_y_continuous(label=scientific_10)+
      scale_x_continuous(label=scientific_10)+
      theme_classic()+theme(legend.position="none"),
    
    legend_plots, nrow=3,ncol=2)
  dev.off()

# ACTIVE DISPERSERS
# Run SpaTemp_BiolComparison.R with the filter for active dispersers
plots_active <- list(plots_HOB_BDD_total[[1]][[3]],
                     plots_HOB_BDD_total[[2]][[3]])

# PASSIVE DISPERSERS
# Run SpaTemp_BiolComparison.R with the filter for passive dispersers
plots_passive <- list(plots_HOB_BDD_total[[1]][[3]],
                      plots_HOB_BDD_total[[4]][[1]],
                      plots_HOB_BDD_total[[4]][[3]])
  
  