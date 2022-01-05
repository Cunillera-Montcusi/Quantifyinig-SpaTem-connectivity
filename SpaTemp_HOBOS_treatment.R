
# Charging Packages ####
library(maps) 
library(mapdata) 
library(maptools)
library(mapproj)
library(rgdal)
library(ggmap)
library(leaflet)
library(tigris)
library(sp)
library(ggplot2)
library(plyr)
library(animation)
library(gridExtra)
library(psych)
library(rstudioapi)
library(data.table)
library(sf)
library(geosphere) 
library(raster)
library(eyelinker)
library(PBSmapping)
library(igraph)
library(adehabitatHR)
library(rgeos)
library(shp2graph)
library(sna)
library(RANN)
library(tidyverse)
library(viridis)

#__________________________________________________________________________
# Charging & depurating HOBOS Dataset ####

setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

Sites <- read.csv("Longlat_Rius.csv", header = T, sep = ";")
colnames(Sites) <- c("Riera", "Codi_HOBO","Latitud","Longitud")

#Correction for matching HOBOS -- collecting the  coordinates of the HOBOS that we have  
for (hob in 1:nrow(Sites)) {
  Matching_HOBOS<- which(as.matrix(dist(Sites[,3:4]))[,hob]==0)
  Random_addition <- rep(seq(0.5,3.5,0.5),nrow(Sites))
  Sites[hob,3:4] <- Sites[hob,3:4]+Random_addition[hob]
}
#Checking
which(as.matrix(dist(Sites[,3:4]))[,1]==0)
plot(Sites$Latitud, Sites$Longitud)

#Charge HOBOS database
HOBOS_sites <- list(
  read.csv("CA_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("M_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("R_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("SA_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("SC_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("T_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("VH_HOBOS_data.csv", header = T, sep = ";"))

# PLotting HOBOS altogether (make the plot window bigger)
par(mfrow=c(3,3))
local_hobos_list <- list()
for (sit in 1:length(HOBOS_sites)) {
  colnames(HOBOS_sites[[sit]])[1] <- c("Day")
  names_hobos <- colnames(HOBOS_sites[[sit]])[2:length(colnames(HOBOS_sites[[sit]]))]
  local_hobos <- c()
  for (nam in 1:length(names_hobos)) {
    local_hobos[nam] <- which(Sites$Codi_HOBO==names_hobos[nam])  
  }
  local_hobos_list[[sit]] <- local_hobos
  plot(Sites$Latitud[local_hobos], Sites$Longitud[local_hobos])
}
par(mfrow=c(1,1))

Sites_list <- list(
  Sites[local_hobos_list[[1]],],
  Sites[local_hobos_list[[2]],],
  Sites[local_hobos_list[[3]],],
  Sites[local_hobos_list[[4]],],
  Sites[local_hobos_list[[5]],],
  Sites[local_hobos_list[[6]],],
  Sites[local_hobos_list[[7]],])

####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
# DIRECTED river network ####
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
# Simple river network ####
####_______________________________________________________________________
# Building river networks based on just directional network.
Simple_river_network <- list()
Simple_river_network_maps <- list()
for (river in 1:length(HOBOS_sites)) {
ST_matrix_out <- matrix(nrow = ncol(HOBOS_sites[[river]])-1,ncol = ncol(HOBOS_sites[[river]])-1, data=0)
spa_connections <-seq(1,ncol(HOBOS_sites[[river]]),1)
time_step_1 <- rep(1,ncol(HOBOS_sites[[river]])-1)
for (site_step in 1:c(ncol(HOBOS_sites[[river]])-2)) {
  #Simple spatial links _______________________
  if(time_step_1[site_step]==1){
    ST_matrix_out[spa_connections[site_step],
                  spa_connections[site_step]+1] <- 1
  }else{
    ST_matrix_out[spa_connections[site_step],
                  spa_connections[site_step]+1] <- 0
    }
  }
Simple_river_network[[river]] <- ST_matrix_out

library(sna)
library(ggnetwork)
#for (e in 1:length(cordenades_xarxes)) {
#  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
#  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
#  CC_values <- PCA_network_results[[e]]

n<- network(Simple_river_network[[river]], directed=T, diag=T)
Simple_river_network_maps[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                              aes(x = x, y = y, xend = xend, yend = yend))+
                                              geom_edges(color = "black", size=1, arrow=arrow(angle = 20)) +
                                              geom_nodes(fill="red", size=5, color="black" ,shape=21)+
                                              theme_classic()+
                                              theme(axis.text = element_blank(),
                                                    axis.ticks = element_blank(),
                                                    legend.position = "none",
                                                    panel.background=element_blank())
                                            
}
png(filename = "Figure/Simple_river_network_maps.png",
    width = 715*6,height = 448*6,units = "px",res = 300)
grid.arrange(Simple_river_network_maps[[1]],Simple_river_network_maps[[2]],Simple_river_network_maps[[3]],
             Simple_river_network_maps[[4]],Simple_river_network_maps[[5]],Simple_river_network_maps[[6]],
             Simple_river_network_maps[[7]], top="Simple river network")
dev.off()


####_______________________________________________________________________
# Simple river network matrix BUILDING ####
####_______________________________________________________________________
direction <- "directed"
weighting <- FALSE
value_LINK <- 1
value_NO_link <- 0

detach("package:sna", unload = TRUE)
# Below there is the function who builds the MATRIX ponderating SPATIAL lINKS=1 and TEMPORAL LINKS=1
ST_matrix_rivers <- list()
ST_directed_Ocloseness_rivers <- list()
ST_directed_Allcloseness_rivers <- list()
ST_directed_betweennes_rivers <- list()


library(doParallel)
registerDoParallel(cores = detectCores())

out_Matrix_LIST <- list()
detach("package:sna", unload = TRUE)
out_Matrix_LIST <- foreach(river=1:length(HOBOS_sites))%dopar%{
#for (river in 1:length(HOBOS_sites)) { - With this it takes 6'26''
# We calculate the number of nodes of our network (used along the function)  
numn_nodes <- ncol(HOBOS_sites[[river]])-1
# We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
### This matrix is the "giant" themplate where we will put all the values.
ST_matrix <- matrix(nrow = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes,
                    ncol = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes, data=0)

ST_Oclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                               numn_nodes, data=0)
ST_Allclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                numn_nodes, data=0)
ST_betweennes_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                               numn_nodes, data=0)

# Once created the template we start to fill it for every day
### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
for (days in 1:(length(HOBOS_sites[[river]]$Day)-1)) {

# First we define the spatial connections of the matrix
### Also known as the rows or columns at which we have to add the values of the connections 
spa_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days-1)*numn_nodes)

# We obtain the time steps:
## time_step_1 is the present
## time_step_2 is the following step (the close future)
time_step_1 <- HOBOS_sites[[river]][days,2:ncol(HOBOS_sites[[river]])]
time_step_2 <- HOBOS_sites[[river]][days+1,2:ncol(HOBOS_sites[[river]])]

#Simple fluvial network_______________________
## This step fills "the diagonal" of each time_step following the direction of the river
## it basically connects the river in a dendritic structure.
for (site_step in 1:c(length(time_step_1)-1)) {
if(time_step_1[site_step]==1){
    ST_matrix[spa_connections[site_step],
              spa_connections[site_step]+1] <- 1 
}else{
    ST_matrix[spa_connections[site_step],
              spa_connections[site_step]+1] <- 0
}}

#FLuvial spatial links _______________________
# Now the party begins. 
## Here we fill the matrix section corresponding to the time_step based on the river graph based on a dendrític. 
require(igraph)
# We create the graph
a <- graph.adjacency(ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                               spa_connections[1]:spa_connections[numn_nodes]], mode=direction,diag = FALSE)

ST_Oclosenness_matrix[days,] <- closeness(a, mode = "out",normalized = T)
ST_Oclosenness_matrix[days,which(time_step_1==0)] <- 0

ST_Allclosenness_matrix[days,] <- closeness(a, mode = "all",normalized = T)

ST_betweennes_matrix[days,] <- betweenness(a)

# We create the matrix where we will drop the information of the shortest paths.
## We will fill "1" or "0" according to the shortest paths. 
All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)
All_river_paths[upper.tri(All_river_paths)] <- value_NO_link

# For each path (e.g., from node 1 to node 7) we "check" the length of the shortest path. 
## check = 0 means that the graph is disconnected.
## check bigger than 0 means that the graph is connected.
for (every_path in 1:c(length(time_step_1)-1)){
  check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
  if (check==0) {
    site <-0
  }else{
  # If bigger than 0. we create a sequence from the path to downstream.
    site <- seq(c(every_path+1),check+every_path,1)
  }

# We fill the "All_river_paths" with 1 on the connections concerning to each "row" or node.
## Site is the vector with the connections (follwing the river downstream).
## When "0" site does not correspond to any row... so the "1" does not go anywhere. 
  All_river_paths[every_path,site] <- value_LINK
}

# We add the "All_river_paths" filled for each node in the "big" matrix specific sites
ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
          spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths

# In the following lines we continue the party towards temporal steps
for (site_step in 1:c(length(time_step_1)-1)) {

#Temporal direct links _______________________
## We generate the temporal connectins
temp_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days)*numn_nodes)  
## We then evaluate the difference between the two time steps and therefore we quantify:
# - Stable links: 0 (WARNING: stable links can be stable 1-1 or 0-0!)
# - Lost links: 1
# - Gained links: -1
## - This values 0,1,-1 define what we will do with the links that match such pattern
temp_change <- time_step_1[site_step]-time_step_2[site_step]

#Stable links (when temp_change=0)
## The mechanics is the same as previously.
if(temp_change==0){# Temporal change is constant 
  if(time_step_1[site_step]==1){# This temporal change implies going from 1 to 1 (so a real stable connected link)
    # We created "All_river_paths"
    All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)
    All_river_paths[upper.tri(All_river_paths)] <- value_NO_link
    # We then do the same as before, check, substitute and add "1" or 0 depending if the connection 
    # following the river is flow.
    for (every_path in site_step:c(length(time_step_1)-1)){
      check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
      if (check==0) {
        site <-0
      }else{
        site <- seq(c(every_path+1),check+every_path,1)
        }
      All_river_paths[every_path,site] <- value_LINK
    }
    # TEMPORAL LINKS are filled in the "future" of our current matrix. This means that we are filling the matrix in 
    # in the diagonal of our "time step" for spatial links but we add the temporal links in the following time step. 
    # so, we evaluate here the present (time step 1) and the future (time step 2) but register as the past (at time step 2)
    ST_matrix[temp_connections[site_step],
              c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
    # Here we add the temporal "link" between "himself". If the link is stable and connected (from 1 to 1), we fill the 
    # diagonal value accordingly. Therefore, we will be able to evaluate the relationship between "himself". Kind of 
    # Tot_Num indicator.  
    ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+numn_nodes] <- value_LINK #1
  }else{# Here we check if the temporal change implies going from 0 to 0 (so a stable disconnected link). Then we put 0
    ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+numn_nodes] <- value_NO_link #0
    }
}

#Lost links (when temp_change=1)
## This just needs to be filled with zeros... so no need to use "All_river_paths"
if(temp_change==1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes] <- value_NO_link #0
}
#Gained links (when temp_change=-1)
## It is a "gain" but it means that "in the present" (time step 1), the node is still disconnected. So it =0
if(temp_change==-1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes]  <- value_NO_link #0
}

#Temporal indirect links _______________________
## Indirect links are those links defined here as the ones that are "lost for the first time" (so temp_change=1).
## When this occurs we asign an "extra" 1 in that particular case. Assuming that when the node dries for the first time
## there is an increase in "dispersal" (downstream directed).
### This only occurs when there is a loss of a previously wet node (from 1 in the present to 0 in the future).
if(temp_change==1){
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
  All_river_paths[upper.tri(All_river_paths)] <- value_NO_link
   for (every_path in site_step:c(length(time_step_1)-1)){
    check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
    if (check==0) {
      site <-0
    }else{
      site <- seq(c(every_path+1),check+every_path,1)
      }
    All_river_paths[every_path,site] <- value_LINK #1
  }
  ST_matrix[temp_connections[site_step],
            c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
  
    }
  }
}

out_Matrix <- list(ST_matrix,ST_Oclosenness_matrix,ST_Allclosenness_matrix,ST_betweennes_matrix)
out_Matrix_LIST[[river]] <- out_Matrix
}


ST_matrix_rivers <- list(out_Matrix_LIST[[1]][[1]],out_Matrix_LIST[[2]][[1]],out_Matrix_LIST[[3]][[1]],
                         out_Matrix_LIST[[4]][[1]],out_Matrix_LIST[[5]][[1]],
                         out_Matrix_LIST[[6]][[1]],out_Matrix_LIST[[7]][[1]])

ST_directed_Ocloseness_rivers <- list(out_Matrix_LIST[[1]][[2]],out_Matrix_LIST[[2]][[2]],out_Matrix_LIST[[3]][[2]],
                                      out_Matrix_LIST[[4]][[2]],out_Matrix_LIST[[5]][[2]],
                                      out_Matrix_LIST[[6]][[2]],out_Matrix_LIST[[7]][[2]])

ST_directed_Allcloseness_rivers <- list(out_Matrix_LIST[[1]][[3]],out_Matrix_LIST[[2]][[3]],out_Matrix_LIST[[3]][[3]],
                                        out_Matrix_LIST[[4]][[3]],out_Matrix_LIST[[5]][[3]],
                                        out_Matrix_LIST[[6]][[3]],out_Matrix_LIST[[7]][[3]])

ST_directed_betweennes_rivers <- list(out_Matrix_LIST[[1]][[4]],out_Matrix_LIST[[2]][[4]],out_Matrix_LIST[[3]][[4]],
                                      out_Matrix_LIST[[4]][[4]],out_Matrix_LIST[[5]][[4]],
                                      out_Matrix_LIST[[6]][[4]],out_Matrix_LIST[[7]][[4]])




####_______________________________________________________________________
# SpatioTemporal matrix collapse calculaiton ####
####_______________________________________________________________________
# Find below the lines to calculate the "collapsing" matrix that just summs all the values of all the SPATIOTEMPORAL matrix 
ST_matrix_plots <- list()
ST_matrix_out_out <- list()

for (river in 1:length(HOBOS_sites)) {
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  # We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
  out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
  for (w in 1:c(length(HOBOS_sites[[river]]$Day)-1)) {
    # We generate the "gen_connections" which corresponds to the columns and rows that we are going to move to obtain and 
    # collapse all the 1, 0 that we have been adding to the main ST matrix
    gen_connections <-seq(1,c(length(colnames(HOBOS_sites[[river]]))-1),1)+w*numn_nodes
    
    # Folloing the "gen_connections" we just "sum" all the present (spatial) and past (temporal) links of the main matrix.  
    out <- ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                     gen_connections[1]:gen_connections[length(gen_connections)]]+
      ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                c(gen_connections[1]+numn_nodes):c(gen_connections[length(gen_connections)]+numn_nodes)]
    
    # Finally we summ all this values in the "common" simple matrix. Here is wher we collapse (sum) all the numbers in one matrix. 
    # The total values should be for a "full" permanent (so always stable links 1 to 1) site the number of total days*2.
    out_out <- out_out+out
    
  }
  
  # We save the collapsed matrix
  ST_matrix_out_out[[river]] <- out_out
  
  # Following lines are just a plotting schema to obtain the graphic representation. 
  # Note that the values are "scaled" to one for colours and sizes.
  n<- network(out_out, directed=T, diag=F)
  
  
  edge_col <- rep(0, length(unlist(n$oel)))
  
  for (edg in 1:c(ncol(HOBOS_sites[[river]])-2)) {
    no_diag_out_out <- out_out
    diag(no_diag_out_out) <- 0
    # This line is a bit tricky but key for coloring properly the edges according to their value! 
    edge_col[n$oel[[edg]]] <- rev(c(no_diag_out_out[edg,]/max(no_diag_out_out))[-seq(edg,1)])
  }
  
  n %e% "Con_values" <- edge_col
  n %e% "Con_values_SIZE" <- ((edge_col)^5)*2
  nod_fill_size <- diag(out_out)/max(diag(out_out))
  n %v% "Site_values" <- nod_fill_size
  ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(colour=Con_values, size=Con_values_SIZE), arrow=arrow(angle = 20), curvature = 0.15) +
    geom_nodes(aes(fill=Site_values, size=Site_values*10), color="black" ,shape=21)+
    scale_color_gradient2(low = "brown",high = "blue",midpoint = 0.5)+
    scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/Dir_NonW_STconnectivityMATRIX.png",
    width = 715*6,height = 448*6,units = "px",res = 300)
grid.arrange(ST_matrix_plots[[1]],ST_matrix_plots[[2]],ST_matrix_plots[[3]],
             ST_matrix_plots[[4]],ST_matrix_plots[[5]],ST_matrix_plots[[6]],
             ST_matrix_plots[[7]], top="ST Dir NonW connectivity Matrix")
dev.off()

####_______________________________________________________________________
# SpatioTemporal connectivity calculaiton ####
####_______________________________________________________________________
ST_connectivity_value <- list()
ST_connectivity_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  # We already know this value
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  # We extract the main matrix that corresponds to the river
  ST_matrix<- ST_matrix_rivers[[river]][1:c(numn_nodes*length(HOBOS_sites[[river]]$Day)),
                                        1:c(numn_nodes*length(HOBOS_sites[[river]]$Day))]
  
  spt_conn <- c()
  # "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
  # considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
  leng_correct <- seq(numn_nodes,1,-1)
  
  # Here we do a several of sum to obtain the indicator: 
  ## First we do an "apply" to sum all the rows. WHich means that we summ all spatial and temporal connections. 
  ## - This number, represents de total amount of times that a node is connected spatially and temporally (max value is num_nodes*2-1)
  ## - We consider that the node can not be connected with itself spatially (this why we have a -1).
  ## Second, we summ all the same nodes together. We summ the results of the previous "apply" to obtain the sum of total amount of times that a node
  ## will be connected to all his possible neighbours through time. 
  ## Third, we correct for the "leng_correct" and make "upstream" values comparable with "downstream" 
  
  out <- foreach(nodes=1:numn_nodes)%dopar%{
  #for (nodes in 1:numn_nodes) {
    spt_conn[nodes] <- sum(apply(ST_matrix,1,sum)[seq(nodes,numn_nodes*length(HOBOS_sites[[river]]$Day)-numn_nodes,numn_nodes)])/leng_correct[nodes]
  }
  # We divide by the number of days so we obtain the "per day" values. max (((numn_nodes*2-1)*512)/position of the node)/513
  spt_conn <- unlist(out)
  spt_conn<- spt_conn/c(length(HOBOS_sites[[river]]$Day))
  
  ST_connectivity_value[[river]] <- spt_conn  
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- c(ST_connectivity_value[[river]][1:(numn_nodes-1)],min(ST_connectivity_value[[river]][1:(numn_nodes-1)])) 
  n %v% "edges_CC_values" <- ST_connectivity_value[[river]][1:(numn_nodes-1)]
  
  ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                          aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/Dir_NonW_STconnectivity.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_connectivity_plot[[1]],ST_connectivity_plot[[2]],ST_connectivity_plot[[3]],
             ST_connectivity_plot[[4]],ST_connectivity_plot[[5]],ST_connectivity_plot[[6]],
             ST_connectivity_plot[[7]], top="ST Dir NonW connectivity")
dev.off()

####_______________________________________________________________________
# SpatioTemporal Out closeness calculaiton ####
####_______________________________________________________________________
ST_Oclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,mean)
  sd_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,sd)
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- mean_Oclo
  n %v% "Sd_values" <- mean_Oclo+sd_Oclo
  n %v% "edges_CC_values" <- mean_Oclo[1:length(mean_Oclo)-1]
  
  ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/Dir_NonW_STOclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_Oclo_plot[[1]],ST_Oclo_plot[[2]],ST_Oclo_plot[[3]],
             ST_Oclo_plot[[4]],ST_Oclo_plot[[5]],ST_Oclo_plot[[6]],
             ST_Oclo_plot[[7]], top="ST Dir NonW Out closennes")
dev.off()

####_______________________________________________________________________
# SpatioTemporal All closeness calculaiton ####
####_______________________________________________________________________
ST_Allclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,mean)
  sd_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,sd)
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- mean_Allclo
  n %v% "Sd_values" <- mean_Allclo+sd_Allclo
  n %v% "edges_CC_values" <- mean_Allclo[1:length(mean_Allclo)-1]
  
  ST_Allclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/Dir_NonW_STAllclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_Allclo_plot[[1]],ST_Allclo_plot[[2]],ST_Allclo_plot[[3]],
             ST_Allclo_plot[[4]],ST_Allclo_plot[[5]],ST_Allclo_plot[[6]],
             ST_Allclo_plot[[7]], top="ST Dir NonW All closennes")
dev.off()

####_______________________________________________________________________
# SpatioTemporal Betweenness calculaiton ####
####_______________________________________________________________________
ST_betw_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_betw <- apply(ST_directed_betweennes_rivers[[river]],2,mean)
  sd_betw <- apply(ST_directed_betweennes_rivers[[river]],2,sd)
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "B_values" <- mean_betw
  n %v% "Sd_values" <- mean_betw+sd_betw
  n %v% "edges_BC_values" <- mean_betw[1:length(mean_betw)-1]
  
  ST_betw_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_BC_values, size=edges_BC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=B_values, size=B_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/Dir_NonW_STbetweenness.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_betw_plot[[1]],ST_betw_plot[[2]],ST_betw_plot[[3]],
             ST_betw_plot[[4]],ST_betw_plot[[5]],ST_betw_plot[[6]],
             ST_betw_plot[[7]], top="ST Dir NonW Betweenness")
dev.off()
####_______________________________________________________________________
# Directed NONWEIGHTED river network OUTPUTS _______________________####

# Global matrix
NonW_ST_matrix_rivers <- ST_matrix_rivers

# Spatiotemporal matrix 
NonW_ST_matrix_out_out <- ST_matrix_out_out
NonW_ST_matrix_plots <- ST_matrix_plots

# Spatiotemporal connectivity 
NonW_ST_connectivity_value <- ST_connectivity_value
NonW_ST_connectivity_plot <- ST_connectivity_plot

# Spatiotemporal Out.closenness 
NonW_ST_directed_Ocloseness_rivers <- ST_directed_Ocloseness_rivers
NonW_ST_Oclo_plot <- ST_Oclo_plot

# Spatiotemporal Out.closenness 
NonW_ST_directed_Allcloseness_rivers <- ST_directed_Allcloseness_rivers
NonW_ST_Allclo_plot <- ST_Allclo_plot

# Spatiotemporal Betweenness 
NonW_ST_directed_betweennes_rivers <- ST_directed_betweennes_rivers
NonW_ST_betw_plot <- ST_betw_plot 

####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________

####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________

# DIRECTED WEIHGTED river network ####
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
#Simple river network weighted matrix BUILDING ####
direction <- "directed"
weighting <- TRUE
value_LINK <- 0
value_NO_link <- 1

detach("package:sna", unload = TRUE)
# Below there is the function who builds the MATRIX ponderating SPATIAL lINKS=1 and TEMPORAL LINKS=1
ST_matrix_rivers <- list()
ST_directed_Ocloseness_rivers <- list()
ST_directed_Allcloseness_rivers <- list()
ST_directed_betweennes_rivers <- list()

#Test final UY_ ALPS 
library(doParallel)
registerDoParallel(cores = detectCores())

out_Matrix_LIST <- list()
detach("package:sna", unload = TRUE)
out_Matrix_LIST <- foreach(river=1:length(HOBOS_sites))%dopar%{
#for (river in 1:length(HOBOS_sites)) {
# We calculate the number of nodes of our network (used along the function)  
numn_nodes <- ncol(HOBOS_sites[[river]])-1
# We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
### This matrix is the "giant" themplate where we will put all the values.
ST_matrix <- matrix(nrow = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes,
                    ncol = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes, data=0)

ST_Oclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                numn_nodes, data=0)
ST_Allclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                  numn_nodes, data=0)
ST_betweennes_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                               numn_nodes, data=0)

# Once created the template we start to fill it for every day
### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
for (days in 1:(length(HOBOS_sites[[river]]$Day)-1)) {
  
  # First we define the spatial connections of the matrix
  ### Also known as the rows or columns at which we have to add the values of the connections 
  spa_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days-1)*numn_nodes)
  
  # We obtain the time steps:
  ## time_step_1 is the present
  ## time_step_2 is the following step (the close future)
  time_step_1 <- HOBOS_sites[[river]][days,2:ncol(HOBOS_sites[[river]])]
  time_step_2 <- HOBOS_sites[[river]][days+1,2:ncol(HOBOS_sites[[river]])]
  
  #Simple fluvial network_______________________
  ## This step fills "the diagonal" of each time_step following the direction of the river
  ## it basically connects the river in a dendritic structure.
if (weighting==T) {
for (site_step in 1:c(length(time_step_1)-1)) {
      if(time_step_1[site_step]==1){
        ST_matrix[spa_connections[site_step],
                  spa_connections[site_step]+1] <- 1*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,site_step+1]
      }else{
        ST_matrix[spa_connections[site_step],
                  spa_connections[site_step]+1] <- 0*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,site_step+1]
      }}  
}else{
  for (site_step in 1:c(length(time_step_1)-1)) {
    if(time_step_1[site_step]==1){
      ST_matrix[spa_connections[site_step],
                spa_connections[site_step]+1] <- value_LINK 
    }else{
      ST_matrix[spa_connections[site_step],
                spa_connections[site_step]+1] <- value_NO_link
    }}  
}
  
#FLuvial spatial links _______________________
# Now the party begins. 
## Here we fill the matrix section corresponding to the time_step based on the river graph based on a dendrític. 
require(igraph)
# We create the graph
a <- graph.adjacency(ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                               spa_connections[1]:spa_connections[numn_nodes]], mode=direction,diag = FALSE,weighted = T)

ST_Oclosenness_matrix[days,] <- closeness(a, mode = "out",normalized = T)
ST_Oclosenness_matrix[days,which(time_step_1==0)] <- 0

ST_Allclosenness_matrix[days,] <- closeness(a, mode = "all",normalized = T)
ST_betweennes_matrix[days,] <- betweenness(a)

# We create the matrix where we will drop the information of the shortest paths.
## We will fill "1" or "0" according to the shortest paths. 
All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)
All_river_paths[upper.tri(All_river_paths)] <- value_NO_link

# For each path (e.g., from node 1 to node 7) we "check" the length of the shortest path. 
## check = 0 means that the graph is disconnected.
## check bigger than 0 means that the graph is connected.
for (every_path in 1:c(length(time_step_1)-1)){
  check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
  if (check==0) {
    site <-0
  }else{
    # If bigger than 0. we create a sequence from the path to downstream.
    site <- seq(c(every_path+1),check+every_path,1)
  }
  
# We fill the "All_river_paths" with 1 on the connections concerning to each "row" or node.
## Site is the vector with the connections (follwing the river downstream).
## When "0" site does not correspond to any row... so the "1" does not go anywhere. 
All_river_paths[every_path,site] <- value_LINK
}
if (weighting==T) {
All_river_paths <- All_river_paths*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))
}
# We add the "All_river_paths" filled for each node in the "big" matrix specific sites
ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
          spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths

# In the following lines we continue the party towards temporal steps
for (site_step in 1:c(length(time_step_1)-1)) {
  
#Temporal direct links _______________________
## We generate the temporal connectins
temp_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days)*numn_nodes)  
## We then evaluate the difference between the two time steps and therefore we quantify:
# - Stable links: 0 (WARNING: stable links can be stable 1-1 or 0-0!)
# - Lost links: 1
# - Gained links: -1
## - This values 0,1,-1 define what we will do with the links that match such pattern
temp_change <- time_step_1[site_step]-time_step_2[site_step]

#Stable links (when temp_change=0)
## The mechanics is the same as previously.
if(temp_change==0){# Temporal change is constant 
  if(time_step_1[site_step]==1){# This temporal change implies going from 1 to 1 (so a real stable connected link)
    # We created "All_river_paths"
    All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)
    All_river_paths[upper.tri(All_river_paths)] <- value_NO_link
    # We then do the same as before, check, substitute and add "1" or 0 depending if the connection 
    # following the river is flow.
    for (every_path in site_step:c(length(time_step_1)-1)){
      check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
      if (check==0) {
      site <-0
    }else{
      site <- seq(c(every_path+1),check+every_path,1)
      }
  All_river_paths[every_path,site] <- value_LINK
}
if (weighting==T) {
All_river_paths <- All_river_paths*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))
}

# TEMPORAL LINKS are filled in the "future" of our current matrix. This means that we are filling the matrix in 
# in the diagonal of our "time step" for spatial links but we add the temporal links in the following time step. 
# so, we evaluate here the present (time step 1) and the future (time step 2) but register as the past (at time step 2)
ST_matrix[temp_connections[site_step],
          c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
# Here we add the temporal "link" between "himself". If the link is stable and connected (from 1 to 1), we fill the 
# diagonal value accordingly. Therefore, we will be able to evaluate the relationship between "himself". Kind of 
# Tot_Num indicator.  
ST_matrix[temp_connections[site_step],
          temp_connections[site_step]+numn_nodes] <- value_LINK

}else{# Here we check if the temporal change implies going from 0 to 0 (so a stable disconnected link). Then we put 0
  ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+numn_nodes] <- value_NO_link
  }
}

#Lost links (when temp_change=1)
## This just needs to be filled with zeros... so no need to use "All_river_paths"
if(temp_change==1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes] <- value_NO_link
}
#Gained links (when temp_change=-1)
## It is a "gain" but it means that "in the present" (time step 1), the node is still disconnected. So it =0
if(temp_change==-1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes]  <- value_NO_link
}

#Temporal indirect links _______________________
## Indirect links are those links defined here as the ones that are "lost for the first time" (so temp_change=1).
## When this occurs we asign an "extra" 1 in that particular case. Assuming that when the node dries for the first time
## there is an increase in "dispersal" (downstream directed).
### This only occurs when there is a loss of a previously wet node (from 1 in the present to 0 in the future).
if(temp_change==1){
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)
  All_river_paths[upper.tri(All_river_paths)] <- value_NO_link
  for (every_path in site_step:c(length(time_step_1)-1)){
    check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
    if (check==0) {
      site <-0
    }else{
      site <- seq(c(every_path+1),check+every_path,1)
    }
    All_river_paths[every_path,site] <- value_LINK
  }
  
if (weighting==T) {
    All_river_paths <- All_river_paths*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))
}
  
ST_matrix[temp_connections[site_step],
        c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
  }
 }
}

out_Matrix <- list(ST_matrix,ST_Oclosenness_matrix,ST_Allclosenness_matrix,ST_betweennes_matrix)
out_Matrix_LIST[[river]] <- out_Matrix
}


ST_matrix_rivers <- list(out_Matrix_LIST[[1]][[1]],out_Matrix_LIST[[2]][[1]],out_Matrix_LIST[[3]][[1]],
                         out_Matrix_LIST[[4]][[1]],out_Matrix_LIST[[5]][[1]],
                         out_Matrix_LIST[[6]][[1]],out_Matrix_LIST[[7]][[1]])

ST_directed_Ocloseness_rivers <- list(out_Matrix_LIST[[1]][[2]],out_Matrix_LIST[[2]][[2]],out_Matrix_LIST[[3]][[2]],
                                      out_Matrix_LIST[[4]][[2]],out_Matrix_LIST[[5]][[2]],
                                      out_Matrix_LIST[[6]][[2]],out_Matrix_LIST[[7]][[2]])

ST_directed_Allcloseness_rivers <- list(out_Matrix_LIST[[1]][[3]],out_Matrix_LIST[[2]][[3]],out_Matrix_LIST[[3]][[3]],
                                        out_Matrix_LIST[[4]][[3]],out_Matrix_LIST[[5]][[3]],
                                        out_Matrix_LIST[[6]][[3]],out_Matrix_LIST[[7]][[3]])

ST_directed_betweennes_rivers <- list(out_Matrix_LIST[[1]][[4]],out_Matrix_LIST[[2]][[4]],out_Matrix_LIST[[3]][[4]],
                                      out_Matrix_LIST[[4]][[4]],out_Matrix_LIST[[5]][[4]],
                                      out_Matrix_LIST[[6]][[4]],out_Matrix_LIST[[7]][[4]])

####_______________________________________________________________________
# SpatioTemporal matrix collapse calculaiton ####
####_______________________________________________________________________
# Find below the lines to calculate the "collapsing" matrix that just summs all the values of all the SPATIOTEMPORAL matrix 
ST_matrix_plots <- list()
ST_matrix_out_out <- list()

for (river in 1:length(HOBOS_sites)) {
numn_nodes <- ncol(HOBOS_sites[[river]])-1

# We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
for (w in 1:c(length(HOBOS_sites[[river]]$Day)-1)) {
# We generate the "gen_connections" which corresponds to the columns and rows that we are going to move to obtain and 
# collapse all the 1, 0 that we have been adding to the main ST matrix
gen_connections <-seq(1,c(length(colnames(HOBOS_sites[[river]]))-1),1)+w*numn_nodes

# Folloing the "gen_connections" we just "sum" all the present (spatial) and past (temporal) links of the main matrix.  
out <- ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                 gen_connections[1]:gen_connections[length(gen_connections)]]+
       ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                 c(gen_connections[1]+numn_nodes):c(gen_connections[length(gen_connections)]+numn_nodes)]

# Finally we summ all this values in the "common" simple matrix. Here is wher we collapse (sum) all the numbers in one matrix. 
# The total values should be for a "full" permanent (so always stable links 1 to 1) site the number of total days*2.
out_out <- out_out+out

}

# We save the collapsed matrix
ST_matrix_out_out[[river]] <- out_out

# Following lines are just a plotting schema to obtain the graphic representation. 
# Note that the values are "scaled" to one for colours and sizes.
n<- network(out_out, directed=T, diag=F)

edge_col <- rep(0, length(unlist(n$oel)))

for (edg in 1:c(ncol(HOBOS_sites[[river]])-2)) {
no_diag_out_out <- out_out
diag(no_diag_out_out) <- 0
# This line is a bit tricky but key for coloring properly the edges according to their value! 
edge_col[n$oel[[edg]]] <- rev(c(no_diag_out_out[edg,]/max(no_diag_out_out))[-seq(edg,1)])
}

n %e% "Con_values" <- edge_col
n %e% "Con_values_SIZE" <- ((edge_col)^5)*2
nod_fill_size <- diag(out_out)/max(diag(out_out))
n %v% "Site_values" <- nod_fill_size
ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                 aes(x = x, y = y, xend = xend, yend = yend))+
                            geom_edges(aes(colour=Con_values, size=Con_values_SIZE), arrow=arrow(angle = 20), curvature = 0.15) +
                            geom_nodes(aes(fill=Site_values, size=Site_values*10), color="black" ,shape=21)+
                            scale_color_gradient2(low = "#9C91CC",high = "#3A341B",midpoint = 0.5)+
                            scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
                            theme_classic()+
                            theme(axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  legend.position = "none",
                                  panel.background=element_blank())
}

png(filename = "Figure/W_Dir_STconnectivityMATRIX.png",
    width = 715*6,height = 448*6,units = "px",res = 300)
grid.arrange(ST_matrix_plots[[1]],ST_matrix_plots[[2]],ST_matrix_plots[[3]],
             ST_matrix_plots[[4]],ST_matrix_plots[[5]],ST_matrix_plots[[6]],
             ST_matrix_plots[[7]], top="ST Dir W connectivity Matrix")
dev.off()

####_______________________________________________________________________
# SpatioTemporal connectivity calculaiton ####
####_______________________________________________________________________
ST_connectivity_value <- list()
ST_connectivity_plot <- list()
for (river in 1:length(HOBOS_sites)) {
# We already know this value
numn_nodes <- ncol(HOBOS_sites[[river]])-1

# We extract the main matrix that corresponds to the river
ST_matrix<- ST_matrix_rivers[[river]][1:c(numn_nodes*length(HOBOS_sites[[river]]$Day)),
                                      1:c(numn_nodes*length(HOBOS_sites[[river]]$Day))]

spt_conn <- c()
# "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
# considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
leng_correct <- seq(numn_nodes,1,-1)

# Here we do a several of sum to obtain the indicator: 
## First we do an "apply" to sum all the rows. WHich means that we summ all spatial and temporal connections. 
## - This number, represents de total amount of times that a node is connected spatially and temporally (max value is num_nodes*2-1)
## - We consider that the node can not be connected with itself spatially (this why we have a -1).
## Second, we summ all the same nodes together. We summ the results of the previous "apply" to obtain the sum of total amount of times that a node
## will be connected to all his possible neighbours through time. 
## Third, we correct for the "leng_correct" and make "upstream" values comparable with "downstream" 

out <- foreach(nodes=1:numn_nodes)%dopar%{
#for (nodes in 1:numn_nodes) {
  spt_conn[nodes] <- sum(apply(ST_matrix,1,sum)[seq(nodes,numn_nodes*length(HOBOS_sites[[river]]$Day)-numn_nodes,numn_nodes)])/leng_correct[nodes]
}
# We divide by the number of days so we obtain the "per day" values. max (((numn_nodes*2-1)*512)/position of the node)/513
spt_conn <- unlist(out)
spt_conn<- spt_conn/c(length(HOBOS_sites[[river]]$Day))

ST_connectivity_value[[river]] <- spt_conn  

n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- c(ST_connectivity_value[[river]][1:(numn_nodes-1)],min(ST_connectivity_value[[river]][1:(numn_nodes-1)])) 
  n %v% "edges_CC_values" <- ST_connectivity_value[[river]][1:(numn_nodes-1)]
  
  ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                          aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/W_Dir_STconnectivity.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_connectivity_plot[[1]],ST_connectivity_plot[[2]],ST_connectivity_plot[[3]],
             ST_connectivity_plot[[4]],ST_connectivity_plot[[5]],ST_connectivity_plot[[6]],
             ST_connectivity_plot[[7]], top="ST Dir W connectivity")
dev.off()

####_______________________________________________________________________
# SpatioTemporal closeness calculaiton ####
####_______________________________________________________________________
ST_Oclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,mean)
  sd_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,sd)
  
n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- c(mean_Oclo[1:c(length(mean_Oclo)-1)],0)
  n %v% "Sd_values" <- mean_Oclo+sd_Oclo
  n %v% "edges_CC_values" <- mean_Oclo[1:c(length(mean_Oclo)-1)]
  
  ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                          aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/W_Dir_ST_Outclosennes_WEIG.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_Oclo_plot[[1]],ST_Oclo_plot[[2]],ST_Oclo_plot[[3]],
             ST_Oclo_plot[[4]],ST_Oclo_plot[[5]],ST_Oclo_plot[[6]],
             ST_Oclo_plot[[7]], top="ST Dir W Out closennes")
dev.off()


####_______________________________________________________________________
# SpatioTemporal All closeness calculaiton ####
####_______________________________________________________________________
ST_Allclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,mean)
  sd_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,sd)
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- mean_Allclo
  n %v% "Sd_values" <- mean_Allclo+sd_Allclo
  n %v% "edges_CC_values" <- mean_Allclo[1:c(length(mean_Allclo)-1)]
  
  ST_Allclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                    aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/W_Dir_STAllclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_Allclo_plot[[1]],ST_Allclo_plot[[2]],ST_Allclo_plot[[3]],
             ST_Allclo_plot[[4]],ST_Allclo_plot[[5]],ST_Allclo_plot[[6]],
             ST_Allclo_plot[[7]], top="ST Dir W All closennes")
dev.off()


####_______________________________________________________________________
# SpatioTemporal Betweenness calculaiton ####
####_______________________________________________________________________
ST_betw_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_betw <- apply(ST_directed_betweennes_rivers[[river]],2,mean)
  sd_betw <- apply(ST_directed_betweennes_rivers[[river]],2,sd)
  
  n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "B_values" <- c(mean_betw[1:c(length(mean_betw)-1)],0)
  n %v% "Sd_values" <- mean_betw+sd_betw
  n %v% "edges_BC_values" <- mean_betw[1:c(length(mean_betw)-1)]
  
  ST_betw_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_BC_values, size=edges_BC_values),arrow=arrow(angle = 20),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=B_values, size=B_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/W_Dir_STbetweenness.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(ST_betw_plot[[1]],ST_betw_plot[[2]],ST_betw_plot[[3]],
             ST_betw_plot[[4]],ST_betw_plot[[5]],ST_betw_plot[[6]],
             ST_betw_plot[[7]], top="ST Dir W Betweenness")
dev.off()


# Directed WEIGHTED river network OUTPUTS _______________________####
# Global matrix
WEIG_ST_matrix_rivers <- ST_matrix_rivers

# Spatiotemporal matrix 
WEIG_ST_matrix_out_out <- ST_matrix_out_out
WEIG_ST_matrix_plots <- ST_matrix_plots

# Spatiotemporal connectivity 
WEIG_ST_connectivity_value <- ST_connectivity_value
WEIG_ST_connectivity_plot <- ST_connectivity_plot

# Spatiotemporal Out.closenness 
WEIG_ST_directed_Ocloseness_rivers <- ST_directed_Ocloseness_rivers
WEIG_ST_Oclo_plot <- ST_Oclo_plot

# Spatiotemporal Out.closenness 
WEIG_ST_directed_Allcloseness_rivers <- ST_directed_Allcloseness_rivers
WEIG_ST_Allclo_plot <- ST_Allclo_plot

# Spatiotemporal Betweenness 
WEIG_ST_directed_betweennes_rivers <- ST_directed_betweennes_rivers
WEIG_ST_betw_plot <- ST_betw_plot 

####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________

# UNDIRECTED river network ####
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
# Complete river network ####
####_______________________________________________________________________
Complete_river_network <- list()
Complete_river_network_maps <- list()
for (river in 1:length(HOBOS_sites)) {
  ST_matrix_out <- matrix(nrow = ncol(HOBOS_sites[[river]])-1,ncol = ncol(HOBOS_sites[[river]])-1, data=0)
  spa_connections <-seq(1,ncol(HOBOS_sites[[river]])-1,1)
  time_step_1 <- rep(1,ncol(HOBOS_sites[[river]])-1)
  for (site_step in 1:c(ncol(HOBOS_sites[[river]])-1)) {
    #Complete spatial links _______________________
    if(time_step_1[site_step]==1){
          ST_matrix_out[spa_connections[site_step],
                        c(1:c(length(spa_connections)))[-site_step]] <- 1
    }else{
      ST_matrix_out[spa_connections[site_step],
                    c(1:c(length(spa_connections)))[-site_step]]<- 0
    }
  }
  Complete_river_network[[river]] <- ST_matrix_out
  
  library(sna)
  library(ggnetwork)
  #for (e in 1:length(cordenades_xarxes)) {
  #  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
  #  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
  #  CC_values <- PCA_network_results[[e]]
  
  n<- network(Complete_river_network[[river]], directed=T, diag=T)
  Complete_river_network_maps[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                               aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "black", size=1, arrow=arrow(angle = 20)) +
    geom_nodes(fill="red", size=5, color="black" ,shape=21)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/Complete_river_network_maps.png",
    width = 715*6,height = 448*6,units = "px",res = 300)
grid.arrange(Complete_river_network_maps[[1]],Complete_river_network_maps[[2]],Complete_river_network_maps[[3]],
             Complete_river_network_maps[[4]],Complete_river_network_maps[[5]],Complete_river_network_maps[[6]],
             Complete_river_network_maps[[7]], top="Complete river network")
dev.off()


####_______________________________________________________________________
# Complete river network matrix BUILDING ####
####_______________________________________________________________________
direction <- "undirected"
value_LINK <- 1
value_NO_link <- 0

Un_ST_matrix_rivers <- list()
Un_ST_directed_Ocloseness_rivers<- list()
Un_ST_directed_Allcloseness_rivers<- list()
Un_ST_directed_betweennes_rivers<- list()

library(doParallel)
registerDoParallel(cores = detectCores())

out_Matrix_LIST <- list()
detach("package:sna", unload = TRUE)
out_Matrix_LIST <- foreach(river=1:length(HOBOS_sites))%dopar%{
#for (river in 1:length(HOBOS_sites)) {
# We calculate the number of nodes of the network   
numn_nodes <- ncol(HOBOS_sites[[river]])-1

# We create the big matrix as before  
Un_ST_matrix <- matrix(nrow = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes,
                      ncol = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes, data=0)

Un_ST_Oclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                numn_nodes, data=0)
Un_ST_Allclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                   numn_nodes, data=0)
Un_ST_betweennes_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                               numn_nodes, data=0)

# We begin to run the following lines for each day
for (days in 1:(length(HOBOS_sites[[river]]$Day)-1)) {
# Spatial connections as in before
spa_connections <-seq(1,length(colnames(HOBOS_sites[[river]]))-1,1)+((days-1)*numn_nodes)
# We analyse the two time steps (present= time_step1 and future= time_step2)
time_step_1 <- HOBOS_sites[[river]][days,2:ncol(HOBOS_sites[[river]])]
time_step_2 <- HOBOS_sites[[river]][days+1,2:ncol(HOBOS_sites[[river]])]

#Complete fluvial network_______________________
for (site_step in 1:length(time_step_1)){
  if(time_step_1[site_step]==1){
    
    Un_ST_matrix[spa_connections[site_step],
                  c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 1
  }else{
    Un_ST_matrix[spa_connections[site_step],
              c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 0
  }
}
#Connections on columns eliminated 
# - UNDIRECTED 
Un_ST_matrix[spa_connections[1:numn_nodes],spa_connections[which(time_step_1==0)]] <- 0

#FLuvial spatial links _______________________
require(igraph)
a <- graph.adjacency(Un_ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                           spa_connections[1]:spa_connections[numn_nodes]],mode=direction,diag = FALSE)

Un_ST_Oclosenness_matrix[days,] <- closeness(a, mode = "out",normalized = T)
Un_ST_Oclosenness_matrix[days,which(time_step_1==0)] <- 0

Un_ST_Allclosenness_matrix[days,] <- closeness(a, mode = "all",normalized = T)

Un_ST_betweennes_matrix[days,] <- betweenness(a)

All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 

for (every_path in 1:c(length(time_step_1))){
  check <- length(all_shortest_paths(a, every_path, c(1:numn_nodes)[-every_path])$res)
  if (check==0) {
    All_river_paths[every_path,] <- value_NO_link
  }else{
    connect_loc <- all_shortest_paths(a, every_path, c(1:numn_nodes)[-every_path])$nrgeo
    site <- which(connect_loc==1)
    All_river_paths[every_path,site] <- value_LINK
  }
}
diag(All_river_paths) <- 0

Un_ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
          spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths

for (site_step in 1:c(length(time_step_1))) {
#Temporal direct links _______________________
temp_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days)*numn_nodes)  
temp_change <- time_step_1[site_step]-time_step_2[site_step]
  
#Stable
if(temp_change==0){
   if(time_step_1[site_step]==1){
    All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0)

     check <- length(all_shortest_paths(a, site_step, c(1:numn_nodes))$res)
     
   if (check==0) {
     All_river_paths[site_step,] <- value_NO_link
    }else{
       connect_loc <- all_shortest_paths(a, site_step, c(1:numn_nodes))$nrgeo
       site <- which(connect_loc==1) 
       All_river_paths[site_step,site] <- value_LINK
  }
  
    diag(All_river_paths) <- 0
    Un_ST_matrix[spa_connections[site_step],
                 c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
   
}else{
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
  diag(All_river_paths) <- 0
  Un_ST_matrix[spa_connections[site_step],
               c(temp_connections[1]):c(temp_connections[numn_nodes])]<- All_river_paths[site_step,]
  }
}

#Lost
if(temp_change==1){
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
  diag(All_river_paths) <- 0
  Un_ST_matrix[spa_connections[site_step],
               c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
}

#Gain
if(temp_change==-1){
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
  diag(All_river_paths) <- 0
  Un_ST_matrix[spa_connections[site_step],
               c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
}

#Temporal indirect links _______________________
#Lost
if(temp_change==1){
  
All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link) 

check <- length(all_shortest_paths(a, site_step, c(1:numn_nodes))$res)
if (check==0) {
  All_river_paths[site_step,] <- value_NO_link
}else{
  connect_loc <- all_shortest_paths(a, site_step, c(1:numn_nodes))$nrgeo
  site <- which(connect_loc==1) 
  All_river_paths[site_step,site] <- value_LINK
}
diag(All_river_paths) <- 0
Un_ST_matrix[spa_connections[site_step],
             c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
}


  }
}

out_Matrix <- list(Un_ST_matrix,Un_ST_Oclosenness_matrix,Un_ST_Allclosenness_matrix,Un_ST_betweennes_matrix)
out_Matrix_LIST[[river]] <- out_Matrix
}


Un_ST_matrix_rivers <- list(out_Matrix_LIST[[1]][[1]],out_Matrix_LIST[[2]][[1]],out_Matrix_LIST[[3]][[1]],
                         out_Matrix_LIST[[4]][[1]],out_Matrix_LIST[[5]][[1]],
                         out_Matrix_LIST[[6]][[1]],out_Matrix_LIST[[7]][[1]])

Un_ST_directed_Ocloseness_rivers <- list(out_Matrix_LIST[[1]][[2]],out_Matrix_LIST[[2]][[2]],out_Matrix_LIST[[3]][[2]],
                                      out_Matrix_LIST[[4]][[2]],out_Matrix_LIST[[5]][[2]],
                                      out_Matrix_LIST[[6]][[2]],out_Matrix_LIST[[7]][[2]])

Un_ST_directed_Allcloseness_rivers <- list(out_Matrix_LIST[[1]][[3]],out_Matrix_LIST[[2]][[3]],out_Matrix_LIST[[3]][[3]],
                                        out_Matrix_LIST[[4]][[3]],out_Matrix_LIST[[5]][[3]],
                                        out_Matrix_LIST[[6]][[3]],out_Matrix_LIST[[7]][[3]])

Un_ST_directed_betweennes_rivers <- list(out_Matrix_LIST[[1]][[4]],out_Matrix_LIST[[2]][[4]],out_Matrix_LIST[[3]][[4]],
                                      out_Matrix_LIST[[4]][[4]],out_Matrix_LIST[[5]][[4]],
                                      out_Matrix_LIST[[6]][[4]],out_Matrix_LIST[[7]][[4]])

####_______________________________________________________________________
# SpatioTemporal matrix calculation ####
####_______________________________________________________________________
Un_ST_matrix_plots <- list()
Un_ST_matrix_out_out <- list()

for (river in 1:length(HOBOS_sites)) {
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
  for (w in 1:c(length(HOBOS_sites[[river]]$Day)-1)) {
    gen_connections <-seq(1,c(length(colnames(HOBOS_sites[[river]]))-1),1)+w*numn_nodes
    
    out <- Un_ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                        gen_connections[1]:gen_connections[length(gen_connections)]]+
           Un_ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                      c(gen_connections[1]+numn_nodes):c(gen_connections[length(gen_connections)]+numn_nodes)]
    out_out <- out_out+out
  }
  
  out_out_out <- matrix(ncol = ncol(out_out), nrow = nrow(out_out),0)
  for (rows in 1:numn_nodes) {
   out_out_out[rows,] <- (out_out[rows,]+out_out[,rows])/2
  }
  out_out_out[lower.tri(out_out_out)] <- 0
  Un_ST_matrix_out_out[[river]] <- out_out_out
  
  # Following lines are just a plotting schema to obtain the graphic representation. 
  # Note that the values are "scaled" to one for colours and sizes.
  n<- network(out_out_out, directed=T, diag=F)
  
  edge_col <- rep(0, length(unlist(n$oel)))
  
  for (edg in 1:c(ncol(HOBOS_sites[[river]])-2)) {
    no_diag_out_out <- out_out_out
    diag(no_diag_out_out) <- 0
    # This line is a bit tricky but key for coloring properly the edges according to their value! 
    edge_col[n$oel[[edg]]] <- rev(c(no_diag_out_out[edg,]/max(no_diag_out_out))[-seq(edg,1)])
  }
  
  n %e% "Con_values" <- edge_col
  n %e% "Con_values_SIZE" <- ((edge_col)^5)*2
  Un_ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(colour=Con_values, size=Con_values_SIZE), arrow=arrow(angle = 20, ends = "both"), curvature = 0.15) +
    geom_nodes(size=3, fill="grey40", color="black" ,shape=21)+
    scale_color_gradient2(low = "brown",high = "blue",midpoint = 0.5)+
    scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/UnD_NonW_STconnectivityMATRIX.png",
    width = 715*10,height = 448*10,units = "px",res = 400)
grid.arrange(Un_ST_matrix_plots[[1]],Un_ST_matrix_plots[[2]],Un_ST_matrix_plots[[3]],
             Un_ST_matrix_plots[[4]],Un_ST_matrix_plots[[5]],Un_ST_matrix_plots[[6]],
             Un_ST_matrix_plots[[7]], top="UnD NonW ST connectivity MATRIX")
dev.off()


####_______________________________________________________________________
# SpatioTemporal connectivity calculation ####
####_______________________________________________________________________
library(viridis)
Un_ST_connectivity_value <- list()
Un_ST_connectivity_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  # We already know this value
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  # We extract the main matrix that corresponds to the river
  Un_ST_matrix<- Un_ST_matrix_rivers[[river]][1:c(numn_nodes*length(HOBOS_sites[[river]]$Day)),
                                              1:c(numn_nodes*length(HOBOS_sites[[river]]$Day))]
  
  spt_conn <- c()

  # Here we do a several of sum to obtain the indicator: 
  ## First we do an "apply" to sum all the rows. WHich means that we summ all spatial and temporal connections. 
  ## - This number, represents de total amount of times that a node is connected spatially and temporally (max value is num_nodes*2-1)
  ## - We consider that the node can not be connected with itself spatially (this why we have a -1).
  ## Second, we summ all the same nodes together. We summ the results of the previous "apply" to obtain the sum of total amount of times that a node
  ## will be connected to all his possible neighbours through time. 
  ## Third, we correct for the "leng_correct" and make "upstream" values comparable with "downstream" 
  
  out <- foreach(nodes=1:numn_nodes)%dopar%{
  #for (nodes in 1:numn_nodes) {
    spt_conn[nodes] <- sum(apply(Un_ST_matrix,1,sum)[seq(nodes,numn_nodes*length(HOBOS_sites[[river]]$Day)-numn_nodes,numn_nodes)])/2
  }
  # We divide by the number of days so we obtain the "per day" values. max (((numn_nodes*2-1)*512)/position of the node)/513
  spt_conn <- unlist(out)
  spt_conn<- spt_conn/c(length(HOBOS_sites[[river]]$Day))
  
  Un_ST_connectivity_value[[river]] <- spt_conn  
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  
  n %v% "edges_CC_values" <- Un_ST_connectivity_value[[river]]/max(Un_ST_connectivity_value[[river]])
  
  Un_ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                          aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both"),curvature = 0.15) +
    geom_nodes(aes(fill =edges_CC_values, size=edges_CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/UnD_NonW_STconnectivity.png",
    width = 715*10,height = 448*10,units = "px",res = 400)
grid.arrange(Un_ST_connectivity_plot[[1]],Un_ST_connectivity_plot[[2]],Un_ST_connectivity_plot[[3]],
             Un_ST_connectivity_plot[[4]],Un_ST_connectivity_plot[[5]],Un_ST_connectivity_plot[[6]],
             Un_ST_connectivity_plot[[7]], top="UnD NonW ST connectivity")
dev.off()


####_______________________________________________________________________
# SpatioTemporal closeness calculaiton ####
####_______________________________________________________________________
Un_ST_Oclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Oclo <- apply(Un_ST_directed_Ocloseness_rivers[[river]],2,mean)
  sd_Oclo <- apply(Un_ST_directed_Ocloseness_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "CC_values" <- mean_Oclo
  n %v% "Sd_values" <- mean_Oclo+sd_Oclo
  n %v% "edges_CC_values" <- mean_Oclo[1:c(length(mean_Oclo)-1)]
  
  Un_ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_NonW_STOclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_Oclo_plot[[1]],Un_ST_Oclo_plot[[2]],Un_ST_Oclo_plot[[3]],
             Un_ST_Oclo_plot[[4]],Un_ST_Oclo_plot[[5]],Un_ST_Oclo_plot[[6]],
             Un_ST_Oclo_plot[[7]], top="UnD NonW ST Out closennes")
dev.off()

####_______________________________________________________________________
# SpatioTemporal All closeness calculaiton ####
####_______________________________________________________________________
Un_ST_Allclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Allclo <- apply(Un_ST_directed_Allcloseness_rivers[[river]],2,mean)
  sd_Allclo <- apply(Un_ST_directed_Allcloseness_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "CC_values" <- mean_Allclo
  n %v% "Sd_values" <- mean_Allclo+sd_Allclo
  n %v% "edges_CC_values" <- mean_Allclo[1:c(length(mean_Allclo)-1)]
  
  Un_ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_NonW_STAllclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_Oclo_plot[[1]],Un_ST_Oclo_plot[[2]],Un_ST_Oclo_plot[[3]],
             Un_ST_Oclo_plot[[4]],Un_ST_Oclo_plot[[5]],Un_ST_Oclo_plot[[6]],
             Un_ST_Oclo_plot[[7]], top="UnD NonW ST All closennes")
dev.off()

####_______________________________________________________________________
# SpatioTemporal Betweenness calculaiton ####
####_______________________________________________________________________
Un_ST_betw_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_betw <- apply(Un_ST_directed_betweennes_rivers[[river]],2,mean)
  sd_betw <- apply(Un_ST_directed_betweennes_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "B_values" <- mean_betw
  n %v% "Sd_values" <- mean_betw+sd_betw
  n %v% "edges_BC_values" <- mean_betw[1:c(length(mean_betw)-1)]
  
  Un_ST_betw_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                  aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_BC_values, size=edges_BC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=B_values, size=B_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_NonW_STbetweenness.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_betw_plot[[1]],Un_ST_betw_plot[[2]],Un_ST_betw_plot[[3]],
             Un_ST_betw_plot[[4]],Un_ST_betw_plot[[5]],Un_ST_betw_plot[[6]],
             Un_ST_betw_plot[[7]], top="UnD NonW ST Betweenness")
dev.off()


# Undirected river network OUTPUTS _______________________####
# River network only

# Global matrix
Un_NonW_ST_matrix_rivers <- Un_ST_matrix_rivers

# Spatiotemporal matrix 
Un_NonW_ST_matrix_out_out <- Un_ST_matrix_out_out
Un_NonW_ST_matrix_plots <- Un_ST_matrix_plots

# Spatiotemporal matrix 
Un_NonW_ST_connectivity_value <- Un_ST_connectivity_value
Un_NonW_ST_connectivity_plot<- Un_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_NonW_ST_directed_Ocloseness_rivers <- Un_ST_directed_Ocloseness_rivers
Un_NonW_ST_Oclo_plot <- Un_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_NonW_ST_directed_Allcloseness_rivers <- Un_ST_directed_Allcloseness_rivers
Un_NonW_ST_Allclo_plot <- Un_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_NonW_ST_directed_betweennes_rivers <- Un_ST_directed_betweennes_rivers
Un_NonW_ST_betw_plot <- Un_ST_betw_plot 

####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________

# UNDIRECTED WEIGHTED river network ####
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________
####_______________________________________________________________________

####_______________________________________________________________________
# Complete river network matrix BUILDING ####
####_______________________________________________________________________
direction <- "undirected"
weighting <- TRUE
value_LINK <- 0.1
value_NO_link <- 1

Un_ST_matrix_rivers <- list()
Un_ST_directed_Ocloseness_rivers<- list()
Un_ST_directed_Allcloseness_rivers<- list()
Un_ST_directed_betweennes_rivers<- list()

#Test final UY_ ALPS 
library(doParallel)
registerDoParallel(cores = detectCores())

out_Matrix_LIST <- list()
detach("package:sna", unload = TRUE)
out_Matrix_LIST <- foreach(river=1:length(HOBOS_sites))%dopar%{
#for (river in i:length(HOBOS_sites)) {
    # We calculate the number of nodes of the network   
    numn_nodes <- ncol(HOBOS_sites[[river]])-1
    
    # We create the big matrix as before  
    Un_ST_matrix <- matrix(nrow = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes,
                           ncol = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes, data=0)
    
    
    Un_ST_Oclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                       numn_nodes, data=0)
    Un_ST_Allclosenness_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                         numn_nodes, data=0)
    Un_ST_betweennes_matrix <- matrix(length(HOBOS_sites[[river]]$Day),
                                      numn_nodes, data=0)
    
    # We begin to run the following lines for each day
    for (days in 1:(length(HOBOS_sites[[river]]$Day)-1)) {
      # Spatial connections as in before
      spa_connections <-seq(1,length(colnames(HOBOS_sites[[river]]))-1,1)+((days-1)*numn_nodes)
      # We analyse the two time steps (present= time_step1 and future= time_step2)
      time_step_1 <- HOBOS_sites[[river]][days,2:ncol(HOBOS_sites[[river]])]
      time_step_2 <- HOBOS_sites[[river]][days+1,2:ncol(HOBOS_sites[[river]])]
      
      #Complete fluvial network_______________________
      if (weighting==T) {
        for (site_step in 1:length(time_step_1)) {
          if(time_step_1[site_step]==1){
            Un_ST_matrix[spa_connections[site_step],
                         c(spa_connections[1]:spa_connections[numn_nodes])] <- 1*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
          }else{
            Un_ST_matrix[spa_connections[site_step],
                         c(spa_connections[1]:spa_connections[numn_nodes])]<- 0*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
          }
        }
      }else{
        for (site_step in 1:length(time_step_1)) {
          if(time_step_1[site_step]==1){
            Un_ST_matrix[spa_connections[site_step],
                         c(spa_connections[1]:spa_connections[numn_nodes])] <- 1
          }else{
            Un_ST_matrix[spa_connections[site_step],
                         c(spa_connections[1]:spa_connections[numn_nodes])]<- 0
          }
        }
      }
      
      #Connections on columns eliminated 
      # - UNDIRECTED 
      Un_ST_matrix[spa_connections[1:numn_nodes],spa_connections[which(time_step_1==0)]] <- 0
      
      #FLuvial spatial links _______________________
      require(igraph)
      a <- graph.adjacency(Un_ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                                        spa_connections[1]:spa_connections[numn_nodes]],mode=direction,diag = FALSE,weighted = T)
      
      Un_ST_Oclosenness_matrix[days,] <- closeness(a, mode = "out",normalized = T)
      Un_ST_Oclosenness_matrix[days,which(time_step_1==0)] <- 0
      Un_ST_Allclosenness_matrix[days,] <- closeness(a, mode="all",normalized = T)
      Un_ST_betweennes_matrix[days,] <- betweenness(a)
      
      All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)

      for (every_path in 1:c(length(time_step_1))){
        check <- length(all_shortest_paths(a, every_path, c(1:numn_nodes)[-every_path])$res)
        if (check==0) {
          All_river_paths[every_path,] <- value_NO_link
        }else{
          connect_loc <- all_shortest_paths(a, every_path, c(1:numn_nodes)[-every_path])$nrgeo
          site <- which(connect_loc==1)
          All_river_paths[every_path,site] <- value_LINK
        }
      }
      diag(All_river_paths) <- 0
      
      if (weighting==T) {
        All_river_paths <- All_river_paths*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))
      }
      
      Un_ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                   spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths
      
      for (site_step in 1:c(length(time_step_1))) {
        #Temporal direct links _______________________
        temp_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days)*numn_nodes)  
        temp_change <- time_step_1[site_step]-time_step_2[site_step]
        
        #Stable
        if(temp_change==0){
          if(time_step_1[site_step]==1){
            All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
            
            check <- length(all_shortest_paths(a, site_step, c(1:numn_nodes)[-site_step])$res)
            if (check==0) {
              All_river_paths[site_step,] <- value_NO_link
            }else{
              connect_loc <- all_shortest_paths(a, site_step, c(1:numn_nodes))$nrgeo
              site <- which(connect_loc==1) 
              All_river_paths[site_step,site] <- value_LINK
            }
            diag(All_river_paths) <- 0
            
            if (weighting==T) {
            All_river_paths[site_step,] <- All_river_paths[site_step,]*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
            }
            
            Un_ST_matrix[spa_connections[site_step],
                         c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
            
          }else{
            All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
            All_river_paths[site_step,] <- All_river_paths[site_step,]*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
            diag(All_river_paths) <- 0
            Un_ST_matrix[spa_connections[site_step],
                       c(temp_connections[1]):c(temp_connections[numn_nodes])]<- All_river_paths[site_step,]
          }
        }
        
        #Lost
        if(temp_change==1){
          All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
          All_river_paths[site_step,] <- All_river_paths[site_step,]*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
          diag(All_river_paths) <- 0
          Un_ST_matrix[spa_connections[site_step],
                     c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
        }
        
        #Gain
        if(temp_change==-1){
          All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link)
          All_river_paths[site_step,] <- All_river_paths[site_step,]*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
          diag(All_river_paths) <- 0
          Un_ST_matrix[spa_connections[site_step],
                     c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
        }
        
        #Temporal indirect links _______________________
        #Lost
        if(temp_change==1){
          
          All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_link) 
          
          check <- length(all_shortest_paths(a, site_step, c(1:numn_nodes))$res)
          if (check==0) {
            All_river_paths[site_step,] <- value_NO_link
          }else{
            connect_loc <- all_shortest_paths(a, site_step, c(1:numn_nodes))$nrgeo
            site <- which(connect_loc==1) 
            All_river_paths[site_step,site] <- value_LINK
          }
          diag(All_river_paths) <- 0
          
          if (weighting==T) {
            All_river_paths[site_step,] <- All_river_paths[site_step,]*as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))[site_step,]
          }
          
          Un_ST_matrix[spa_connections[site_step],
                     c(temp_connections[1]):c(temp_connections[numn_nodes])] <- All_river_paths[site_step,]
        }
      }
    }

    out_Matrix <- list(Un_ST_matrix,Un_ST_Oclosenness_matrix,Un_ST_Allclosenness_matrix,Un_ST_betweennes_matrix)
    out_Matrix_LIST[[river]] <- out_Matrix
}


Un_ST_matrix_rivers <- list(out_Matrix_LIST[[1]][[1]],out_Matrix_LIST[[2]][[1]],out_Matrix_LIST[[3]][[1]],
                            out_Matrix_LIST[[4]][[1]],out_Matrix_LIST[[5]][[1]],
                            out_Matrix_LIST[[6]][[1]],out_Matrix_LIST[[7]][[1]])

Un_ST_directed_Ocloseness_rivers <- list(out_Matrix_LIST[[1]][[2]],out_Matrix_LIST[[2]][[2]],out_Matrix_LIST[[3]][[2]],
                                         out_Matrix_LIST[[4]][[2]],out_Matrix_LIST[[5]][[2]],
                                         out_Matrix_LIST[[6]][[2]],out_Matrix_LIST[[7]][[2]])

Un_ST_directed_Allcloseness_rivers <- list(out_Matrix_LIST[[1]][[3]],out_Matrix_LIST[[2]][[3]],out_Matrix_LIST[[3]][[3]],
                                           out_Matrix_LIST[[4]][[3]],out_Matrix_LIST[[5]][[3]],
                                           out_Matrix_LIST[[6]][[3]],out_Matrix_LIST[[7]][[3]])

out_Matrix_LIST[[1]][[4]][1:513,1:15] <- 0

Un_ST_directed_betweennes_rivers <- list(out_Matrix_LIST[[1]][[4]],out_Matrix_LIST[[2]][[4]],out_Matrix_LIST[[3]][[4]],
                                         out_Matrix_LIST[[4]][[4]],out_Matrix_LIST[[5]][[4]],
                                         out_Matrix_LIST[[6]][[4]],out_Matrix_LIST[[7]][[4]])


####_______________________________________________________________________
# SpatioTemporal matrix calculation ####
####_______________________________________________________________________
Un_ST_matrix_plots <- list()
Un_ST_matrix_out_out <- list()

for (river in 1:length(HOBOS_sites)) {
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
  for (w in 1:c(length(HOBOS_sites[[river]]$Day)-1)) {
    gen_connections <-seq(1,c(length(colnames(HOBOS_sites[[river]]))-1),1)+w*numn_nodes
    
    out <- Un_ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                        gen_connections[1]:gen_connections[length(gen_connections)]]+
      Un_ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                                   c(gen_connections[1]+numn_nodes):c(gen_connections[length(gen_connections)]+numn_nodes)]
    out_out <- out_out+out
  }
  
  out_out_out <- matrix(ncol = ncol(out_out), nrow = nrow(out_out),0)
  for (rows in 1:numn_nodes) {
    out_out_out[rows,] <- (out_out[rows,]+out_out[,rows])/2
  }
  out_out_out[lower.tri(out_out_out)] <- 0
  Un_ST_matrix_out_out[[river]] <- out_out_out
  
  # Following lines are just a plotting schema to obtain the graphic representation. 
  # Note that the values are "scaled" to one for colours and sizes.
  n<- network(out_out_out, directed=T, diag=F)
  
  edge_col <- rep(0, length(unlist(n$oel)))
  
  for (edg in 1:c(ncol(HOBOS_sites[[river]])-2)) {
    no_diag_out_out <- out_out_out
    diag(no_diag_out_out) <- 0
    # This line is a bit tricky but key for coloring properly the edges according to their value! 
    edge_col[n$oel[[edg]]] <- rev(c(no_diag_out_out[edg,]/max(no_diag_out_out))[-seq(edg,1)])
  }
  
  n %e% "Con_values" <- edge_col
  n %e% "Con_values_SIZE" <- ((edge_col)^5)*2
  Un_ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                        aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(colour=Con_values, size=Con_values_SIZE), arrow=arrow(angle = 20, ends = "both"), curvature = 0.15) +
    geom_nodes(size=3, fill="grey40", color="black" ,shape=21)+
    scale_color_gradient2(low = "#9C91CC",high = "#3A341B",midpoint = 0.5)+
    scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/UnD_Weight_STconnectivityMATRIX.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_matrix_plots[[1]],Un_ST_matrix_plots[[2]],Un_ST_matrix_plots[[3]],
             Un_ST_matrix_plots[[4]],Un_ST_matrix_plots[[5]],Un_ST_matrix_plots[[6]],
             Un_ST_matrix_plots[[7]], top="ST UnD W connectivity Matrix" )
dev.off()


####_______________________________________________________________________
# SpatioTemporal connectivity calculation ####
####_______________________________________________________________________
library(viridis)
Un_ST_connectivity_value <- list()
Un_ST_connectivity_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  # We already know this value
  numn_nodes <- ncol(HOBOS_sites[[river]])-1
  
  # We extract the main matrix that corresponds to the river
  Un_ST_matrix<- Un_ST_matrix_rivers[[river]][1:c(numn_nodes*length(HOBOS_sites[[river]]$Day)),
                                              1:c(numn_nodes*length(HOBOS_sites[[river]]$Day))]
  
  spt_conn <- c()
  
  # Here we do a several of sum to obtain the indicator: 
  ## First we do an "apply" to sum all the rows. WHich means that we summ all spatial and temporal connections. 
  ## - This number, represents de total amount of times that a node is connected spatially and temporally (max value is num_nodes*2-1)
  ## - We consider that the node can not be connected with itself spatially (this why we have a -1).
  ## Second, we summ all the same nodes together. We summ the results of the previous "apply" to obtain the sum of total amount of times that a node
  ## will be connected to all his possible neighbours through time. 
  ## Third, we correct for the "leng_correct" and make "upstream" values comparable with "downstream" 
  
  out <- foreach(nodes=1:numn_nodes)%dopar%{
  #for (nodes in 1:numn_nodes) {
    spt_conn[nodes] <- sum(apply(Un_ST_matrix,1,sum)[seq(nodes,numn_nodes*length(HOBOS_sites[[river]]$Day)-numn_nodes,numn_nodes)])/2
  }
  # We divide by the number of days so we obtain the "per day" values. max (((numn_nodes*2-1)*512)/position of the node)/513
  spt_conn <- unlist(out)
  spt_conn<- spt_conn/c(length(HOBOS_sites[[river]]$Day))
  
  Un_ST_connectivity_value[[river]] <- spt_conn  
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  
  n %v% "CC_values" <- Un_ST_connectivity_value[[river]] 
  n %v% "edges_CC_values" <- Un_ST_connectivity_value[[river]]/max(Un_ST_connectivity_value[[river]])
  
  Un_ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                             aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both"),curvature = 0.15) +
    geom_nodes(size=3, fill="grey40", color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

png(filename = "Figure/UnD_Weight_STconnectivity.png",
    width = 715*12,height = 448*12,units = "px",res = 400)
grid.arrange(Un_ST_connectivity_plot[[1]],Un_ST_connectivity_plot[[2]],Un_ST_connectivity_plot[[3]],
             Un_ST_connectivity_plot[[4]],Un_ST_connectivity_plot[[5]],Un_ST_connectivity_plot[[6]],
             Un_ST_connectivity_plot[[7]], top="ST UnD W connectivity")
dev.off()


####_______________________________________________________________________
# SpatioTemporal closeness calculaiton ####
####_______________________________________________________________________
Un_ST_Oclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Oclo <- apply(Un_ST_directed_Ocloseness_rivers[[river]],2,mean)
  sd_Oclo <- apply(Un_ST_directed_Ocloseness_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "CC_values" <- mean_Oclo
  n %v% "Sd_values" <- mean_Oclo+sd_Oclo
  n %v% "edges_CC_values" <- mean_Oclo
  
  Un_ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_Weight_STOclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_Oclo_plot[[1]],Un_ST_Oclo_plot[[2]],Un_ST_Oclo_plot[[3]],
             Un_ST_Oclo_plot[[4]],Un_ST_Oclo_plot[[5]],Un_ST_Oclo_plot[[6]],
             Un_ST_Oclo_plot[[7]], top="ST UnD W Out closennes")
dev.off()

####_______________________________________________________________________
# SpatioTemporal All closeness calculaiton ####
####_______________________________________________________________________
Un_ST_Allclo_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_Allclo <- apply(Un_ST_directed_Allcloseness_rivers[[river]],2,mean)
  sd_Allclo <- apply(Un_ST_directed_Allcloseness_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "CC_values" <- mean_Allclo
  n %v% "Sd_values" <- mean_Allclo+sd_Allclo
  n %v% "edges_CC_values" <- mean_Allclo
  
  Un_ST_Allclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_Weight_STAllclosennes.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_Allclo_plot[[1]],Un_ST_Allclo_plot[[2]],Un_ST_Allclo_plot[[3]],
             Un_ST_Allclo_plot[[4]],Un_ST_Allclo_plot[[5]],Un_ST_Allclo_plot[[6]],
             Un_ST_Allclo_plot[[7]], top="ST UnD W All closennes")
dev.off()


####_______________________________________________________________________
# SpatioTemporal Betweenness calculaiton ####
####_______________________________________________________________________
Un_ST_betw_plot <- list()
for (river in 1:length(HOBOS_sites)) {
  mean_betw <- apply(Un_ST_directed_betweennes_rivers[[river]],2,mean)
  sd_betw <- apply(Un_ST_directed_betweennes_rivers[[river]],2,sd)
  
  Network_representation <- Complete_river_network[[river]]
  Network_representation[lower.tri(Network_representation)] <- 0
  n<- network(Network_representation, directed=T, diag=T)
  n %v% "B_values" <- mean_betw
  n %v% "Sd_values" <- mean_betw+sd_betw
  n %v% "edges_BC_values" <- mean_betw
  
  Un_ST_betw_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                     aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_BC_values, size=edges_BC_values/2),arrow=arrow(angle = 20,ends = "both" ),curvature = 0.15) +
    geom_nodes(aes(size=Sd_values),fill="grey50", color="black" ,shape=21, alpha=0.2)+
    geom_nodes(aes(fill=B_values, size=B_values), color="black" ,shape=21)+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
  
}

png(filename = "Figure/UnD_Weight_STbetweenness.png",
    width = 715*6,height = 448*6,units = "px",res = 400)
grid.arrange(Un_ST_betw_plot[[1]],Un_ST_betw_plot[[2]],Un_ST_betw_plot[[3]],
             Un_ST_betw_plot[[4]],Un_ST_betw_plot[[5]],Un_ST_betw_plot[[6]],
             Un_ST_betw_plot[[7]], top="ST UnD W Betweenness")
dev.off()

# Undirected WEIGHTED river network OUTPUTS _______________________####

# Global matrix
Un_WEIG_ST_matrix_rivers <- Un_ST_matrix

# Spatiotemporal matrix 
Un_WEIG_ST_matrix_out_out <- Un_ST_matrix_out_out
Un_WEIG_ST_matrix_plots <- Un_ST_matrix_plots

# Spatiotemporal matrix 
Un_WEIG_ST_connectivity_value <- Un_ST_connectivity_value
Un_WEIG_ST_connectivity_plot<- Un_ST_connectivity_plot

# Spatiotemporal All closenness 
Un_WEIG_ST_directed_Ocloseness_rivers <- Un_ST_directed_Ocloseness_rivers
Un_WEIG_ST_Oclo_plot <- Un_ST_Oclo_plot

# Spatiotemporal Out closenness 
Un_WEIG_ST_directed_Allcloseness_rivers <- Un_ST_directed_Allcloseness_rivers
Un_WEIG_ST_Allclo_plot <- Un_ST_Allclo_plot

# Spatiotemporal Betweenness 
Un_WEIG_ST_directed_betweennes_rivers <- Un_ST_directed_betweennes_rivers
Un_WEIG_ST_betw_plot <- Un_ST_betw_plot 


