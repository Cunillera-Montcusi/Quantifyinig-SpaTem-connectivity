
# Directed NONWEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
NonW_ST_matrix_rivers

# Spatiotemporal matrix 
NonW_ST_matrix_out_out
NonW_ST_directed_matrix_mean <- list()
NonW_ST_directed_matrix_sd <- list()
for (river in 1:length(NonW_ST_matrix_out_out)) {
  NonW_ST_directed_matrix_mean[[river]] <- apply(NonW_ST_matrix_out_out[[river]],1,mean)
  NonW_ST_directed_matrix_sd[[river]] <- apply(NonW_ST_matrix_out_out[[river]],1,sd)
}

# Spatiotemporal connectivity 
NonW_ST_connectivity_value
NonW_ST_connectivity_plot

# Spatiotemporal Out.closenness 
NonW_ST_directed_Ocloseness_rivers
NonW_ST_directed_Oclo_mean <- list()
NonW_ST_directed_Oclo_sd <- list()
for (river in 1:length(NonW_ST_directed_Ocloseness_rivers)) {
  NonW_ST_directed_Oclo_mean[[river]] <- apply(NonW_ST_directed_Ocloseness_rivers[[river]],2,mean)
  NonW_ST_directed_Oclo_sd[[river]] <- apply(NonW_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

NonW_ST_oclo_plot

# Spatiotemporal ALL.closenness 
NonW_ST_directed_Allcloseness_rivers
NonW_ST_directed_Allclo_mean <- list()
NonW_ST_directed_Allclo_sd <- list()
for (river in 1:length(NonW_ST_directed_Allcloseness_rivers)) {
  NonW_ST_directed_Allclo_mean[[river]] <- apply(NonW_ST_directed_Allcloseness_rivers[[river]],2,mean)
  NonW_ST_directed_Allclo_sd[[river]] <- apply(NonW_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

NonW_ST_Allclo_plot

# Spatiotemporal Betweenness 
NonW_ST_directed_betweennes_rivers
NonW_ST_directed_betw_mean <- list()
NonW_ST_directed_betw_sd <- list()
for (river in 1:length(NonW_ST_directed_betweennes_rivers)) {
  NonW_ST_directed_betw_mean[[river]] <- apply(NonW_ST_directed_betweennes_rivers[[river]],2,mean)
  NonW_ST_directed_betw_sd[[river]] <- apply(NonW_ST_directed_betweennes_rivers[[river]],2,sd)
}

NonW_ST_betw_plot


NonW_ST_directed_output <- data.frame(NonW_Dir_mat=unlist(NonW_ST_directed_matrix_mean),
                                      NonW_Dir_con=unlist(NonW_ST_connectivity_value),
                                      NonW_Dir_Ocl=unlist(NonW_ST_directed_Oclo_mean),
                                      NonW_Dir_Acl=unlist(NonW_ST_directed_Allclo_mean),
                                      NonW_Dir_Bet=unlist(NonW_ST_directed_betw_mean))

corrmorant::corrmorant(cbind())

# Directed WEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
WEIG_ST_matrix_rivers 

# Spatiotemporal matrix 
WEIG_ST_matrix_out_out
WEIG_ST_matrix_mean <- list()
WEIG_ST_matrix_sd <- list()
for (river in 1:length(WEIG_ST_matrix_out_out)) {
  WEIG_ST_matrix_mean[[river]] <- apply(WEIG_ST_matrix_out_out[[river]],1,mean)
  WEIG_ST_matrix_sd[[river]] <- apply(WEIG_ST_matrix_out_out[[river]],1,sd)
}

WEIG_ST_matrix_plots

# Spatiotemporal connectivity 
WEIG_ST_connectivity_value
WEIG_ST_connectivity_plot

# Spatiotemporal Out.closenness 
WEIG_ST_directed_Ocloseness_rivers
WEIG_ST_Oclo_mean <- list()
WEIG_ST_Oclo_sd <- list()
for (river in 1:length(WEIG_ST_directed_Ocloseness_rivers)) {
  WEIG_ST_Oclo_mean[[river]] <- apply(WEIG_ST_directed_Ocloseness_rivers[[river]],2,mean)
  WEIG_ST_Oclo_sd[[river]] <- apply(WEIG_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

WEIG_ST_Oclo_plot

# Spatiotemporal All.closenness 
WEIG_ST_directed_Allcloseness_rivers
WEIG_ST_Allclo_mean <- list()
WEIG_ST_Allclo_sd <- list()
for (river in 1:length(WEIG_ST_directed_Allcloseness_rivers)) {
  WEIG_ST_Allclo_mean[[river]] <- apply(WEIG_ST_directed_Allcloseness_rivers[[river]],2,mean)
  WEIG_ST_Allclo_sd[[river]] <- apply(WEIG_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

WEIG_ST_Allclo_plot

# Spatiotemporal Betweenness 
WEIG_ST_directed_betweennes_rivers
WEIG_ST_betw_mean <- list()
WEIG_ST_betw_sd <- list()
for (river in 1:length(WEIG_ST_directed_betweennes_rivers)) {
  WEIG_ST_betw_mean[[river]] <- apply(WEIG_ST_directed_betweennes_rivers[[river]],2,mean)
  WEIG_ST_betw_sd[[river]] <- apply(WEIG_ST_directed_betweennes_rivers[[river]],2,sd)
}

WEIG_ST_betw_plot


WEIG_ST_directed_output <- data.frame(WEIG_Dir_mat=unlist(WEIG_ST_matrix_mean),
                                      WEIG_Dir_con=unlist(WEIG_ST_connectivity_value),
                                      WEIG_Dir_Ocl=unlist(WEIG_ST_Oclo_mean),
                                      WEIG_Dir_Acl=unlist(WEIG_ST_Allclo_mean),
                                      WEIG_Dir_Bet=unlist(WEIG_ST_betw_mean))


# Undirected river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
Un_NonW_ST_matrix_rivers

# Spatiotemporal matrix 
Un_NonW_ST_matrix_out_out
Un_NonW_ST_matrix_mean <- list()
Un_NonW_ST_matrix_sd <- list()
for (river in 1:length(Un_NonW_ST_matrix_out_out)) {
  Un_NonW_ST_matrix_mean[[river]] <- apply(Un_NonW_ST_matrix_out_out[[river]],1,mean)
  Un_NonW_ST_matrix_sd[[river]] <- apply(Un_NonW_ST_matrix_out_out[[river]],1,sd)
}

Un_NonW_ST_matrix_plots

# Spatiotemporal matrix 
Un_NonW_ST_connectivity_value
Un_NonW_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_NonW_ST_directed_Ocloseness_rivers
Un_NonW_ST_Oclo_mean <- list()
Un_NonW_ST_Oclo_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_Ocloseness_rivers)) {
  Un_NonW_ST_Oclo_mean[[river]] <- apply(Un_NonW_ST_directed_Ocloseness_rivers[[river]],2,mean)
  Un_NonW_ST_Oclo_sd[[river]] <- apply(Un_NonW_ST_directed_Ocloseness_rivers[[river]],2,sd)
}

Un_NonW_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_NonW_ST_directed_Allcloseness_rivers
Un_NonW_ST_Allclo_mean <- list()
Un_NonW_ST_Allclo_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_Allcloseness_rivers)) {
  Un_NonW_ST_Allclo_mean[[river]] <- apply(Un_NonW_ST_directed_Allcloseness_rivers[[river]],2,mean)
  Un_NonW_ST_Allclo_sd[[river]] <- apply(Un_NonW_ST_directed_Allcloseness_rivers[[river]],2,sd)
}

Un_NonW_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_NonW_ST_directed_betweennes_rivers
Un_NonW_ST_betw_mean <- list()
Un_NonW_ST_betw_sd <- list()
for (river in 1:length(Un_NonW_ST_directed_betweennes_rivers)) {
  Un_NonW_ST_betw_mean[[river]] <- apply(Un_NonW_ST_directed_betweennes_rivers[[river]],2,mean)
  Un_NonW_ST_betw_sd[[river]] <- apply(Un_NonW_ST_directed_betweennes_rivers[[river]],2,sd)
}
Un_NonW_ST_betw_plot

Un_NonW_ST_output <- data.frame(Un_NonW_mat=unlist(Un_NonW_ST_matrix_mean),
                                Un_NonW_con=unlist(Un_NonW_ST_connectivity_value),
                                Un_NonW_Ocl=unlist(Un_NonW_ST_Oclo_mean),
                                Un_NonW_Acl=unlist(Un_NonW_ST_Allclo_mean),
                                Un_NonW_Bet=unlist(Un_NonW_ST_betw_mean))


# Undirected WEIGHTED river network OUTPUTS ___________________________________________________________________________________________________________####
# Global matrix
Un_WEIG_ST_matrix_rivers <- Un_ST_matrix

# Spatiotemporal matrix 
Un_WEIG_ST_matrix_out_out <- Un_ST_matrix_out_out
Un_WEIG_ST_matrix_mean <- list()
Un_WEIG_ST_matrix_sd <- list()
for (river in 1:length(Un_WEIG_ST_matrix_out_out)) {
  Un_WEIG_ST_matrix_mean[[river]] <- apply(Un_WEIG_ST_matrix_out_out[[river]],1,mean)
  Un_WEIG_ST_matrix_sd[[river]] <- apply(Un_WEIG_ST_matrix_out_out[[river]],1,sd)
}

Un_WEIG_ST_matrix_plots <- Un_ST_matrix_plots

# Spatiotemporal matrix 
Un_WEIG_ST_connectivity_value <- Un_ST_connectivity_value
Un_WEIG_ST_connectivity_plot<- Un_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_WEIG_ST_directed_Ocloseness_rivers <- Un_ST_directed_Ocloseness_rivers
Un_WEIG_ST_Oclo_mean <- list()
Un_WEIG_ST_Oclo_sd <- list()
for (river in 1:length(Un_WEIG_ST_directed_Ocloseness_rivers)) {
  Un_WEIG_ST_Oclo_mean[[river]] <- apply(Un_WEIG_ST_directed_Ocloseness_rivers[[river]],2,mean)
  Un_WEIG_ST_Oclo_sd[[river]] <- apply(Un_WEIG_ST_directed_Ocloseness_rivers[[river]],2,sd)
}
Un_WEIG_ST_Oclo_plot <- Un_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_WEIG_ST_directed_Allcloseness_rivers <- Un_ST_directed_Allcloseness_rivers
Un_WEIG_ST_Allclo_plot <- Un_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_WEIG_ST_directed_betweennes_rivers <- Un_ST_directed_betweennes_rivers
Un_WEIG_ST_betw_plot <- Un_ST_betw_plot 

