
grid.arrange(
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$Dir_NonW_ST_con_plo[[7]],bottom="1.1) Directed Binary"),
  arrangeGrob(Dir_WEIG_Net$Dir_WEIG_ST_con_plo[[7]],bottom="1.2) Directed Weighted"),
  arrangeGrob(UnD_NonW_Net$UnD_NonW_ST_con_plo[[7]],bottom="2.1) Undirected Binary"),
  arrangeGrob(UnD_WEIG_Net$UnD_WEIG_ST_con_plo[[7]],bottom="2.2) Undirected Weighted"),
  top="A) STcon"
),
arrangeGrob(
  arrangeGrob(Dir_NonW_Net$Dir_NonW_ST_matrix_plo[[7]],bottom="1.1) Directed Binary"),
  arrangeGrob(Dir_WEIG_Net$Dir_WEIG_ST_matrix_plo[[7]],bottom="1.2) Directed Weighted"),
  arrangeGrob(UnD_NonW_Net$UnD_NonW_ST_matrix_plo[[7]],bottom="2.1) Undirected Binary"),
  arrangeGrob(UnD_WEIG_Net$UnD_WEIG_ST_matrix_plo[[7]],bottom="2.2) Undirected Weighted"),
  top="B) STconmat"
)
)





data_PCA <- data.frame(data.frame("STcon 1.1."=NonW_ST_directed_out$NonW_Dir_con,
                                  "STcon 1.2."=WEIG_ST_directed_out$WEIG_Dir_con,
                                  "STcon 2.1."=Un_NonW_ST_out$Un_NonW_con, 
                                  "STcon 2.2."=Un_WEIG_ST_out$Un_WEIG_con),
            Old_HOBOS_comb[,3:5])

png(filename =paste("Figure/Data_treat/","PCAs_Fig3",".png"), 
    width = 700*2, height = 600*2, 
    units = "px",res = 300)
autoplot(prcomp(data_PCA,
  center = T, scale. = T),
  loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
  loadings.colour = 'red', loadings.label.colour="black")+
  geom_point(size=2, alpha=0.7,shape=21,
             color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
  scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
  theme_bw()+labs(title =names_approxim[1])
dev.off()

# Printing image for a PCA with all the variables together 
for (a in 1:length(four_approxim)) {
  png(filename =paste("Figure/Data_treat/","PCAs ST ind. data vs TOT",".png"), 
      width = 700*4, height = 600*4, 
      units = "px",res = 300)
  grid.arrange(
    autoplot(prcomp(
      data.frame(NonW_ST_directed_out[,3],Old_HOBOS_comb[,3:5]),
      center = T, scale. = T),
      loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
      loadings.colour = 'red', loadings.label.colour="black")+
      geom_point(size=2, alpha=0.7,shape=21,
                 color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
      scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
      theme_bw()+labs(title =names_approxim[1]),
    
    autoplot(prcomp(
      data.frame(WEIG_ST_directed_out[,3],Old_HOBOS_comb[,3:5]),
      center = T, scale. = T),
      loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
      loadings.colour = 'red', loadings.label.colour="black")+
      geom_point(size=2, alpha=0.7,shape=21,
                 color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
      scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
      theme_bw()+labs(title =names_approxim[2]),
    
    autoplot(prcomp(
      data.frame(Un_NonW_ST_out[,3],Old_HOBOS_comb[,3:5]),
      center = T, scale. = T),
      loadings=T, loadings.label = TRUE, loadings.label.size = 5, shape = FALSE, label=F,
      loadings.colour = 'red', loadings.label.colour="black")+
      geom_point(size=2, alpha=0.7,shape=21,
                 color="grey30", aes(fill=NonW_ST_directed_out[,1]))+
      scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
      theme_bw()+labs(title =names_approxim[3]),
    
    autoplot(prcomp(
      data.frame(Un_WEIG_ST_out[,3],Old_HOBOS_comb[,3:5]),
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
