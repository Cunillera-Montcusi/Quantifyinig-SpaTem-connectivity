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

setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

Sites <- read.csv("Longlat_Rius.csv", header = T, sep = ";")
colnames(Sites) <- c("Riera", "Codi_HOBO","Latitud","Longitud")

#Correction for matching HOBOS
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

grid.arrange(Simple_river_network_maps[[1]],Simple_river_network_maps[[2]],Simple_river_network_maps[[3]],
             Simple_river_network_maps[[4]],Simple_river_network_maps[[5]],Simple_river_network_maps[[6]],
             Simple_river_network_maps[[7]])


# SPATIOREMPORAL CONNECTIVITY NETWORKS

defined_neighbours <- 2

ST_matrix_rivers <- list()
for (river in 1:length(HOBOS_sites)) {
  
numn_nodes <- ncol(HOBOS_sites[[river]])-1
ST_matrix <- matrix(nrow = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes,
                    ncol = numn_nodes*length(HOBOS_sites[[river]]$Day)+numn_nodes, data=0)
for (days in 1:(length(HOBOS_sites[[river]]$Day)-1)) {

spa_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days-1)*numn_nodes)

time_step_1 <- HOBOS_sites[[river]][days,2:ncol(HOBOS_sites[[river]])]
time_step_2 <- HOBOS_sites[[river]][days+1,2:ncol(HOBOS_sites[[river]])]

for (site_step in 1:c(length(time_step_1)-1)) {
#Simple spatial links _______________________
if(time_step_1[site_step]==1){
    ST_matrix[spa_connections[site_step],
              spa_connections[site_step]+1] <- 1
}else{
    ST_matrix[spa_connections[site_step],
              spa_connections[site_step]+1] <- 0
}}
#FLuvial spatial links _______________________
require(igraph)
a <- graph.adjacency(ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                               spa_connections[1]:spa_connections[numn_nodes]],mode="directed",diag = FALSE)
All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
for (every_path in 1:c(length(time_step_1)-1)){
  check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
  if (check==0) {
    site <-0
  }else{
    site <- seq(c(every_path+1),check+every_path,1)
  }
  All_river_paths[every_path,site[1]:site[length(site)]] <- 1
}
ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
          spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths

for (site_step in 1:c(length(time_step_1)-1)) {
#Temporal direct links _______________________
temp_connections <-seq(1,length(colnames(HOBOS_sites[[river]])),1)+((days)*numn_nodes)  

temp_change <- time_step_1[site_step]-time_step_2[site_step]
#Stable
if(temp_change==0){
  if(time_step_1[site_step]==1){

    
    All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
    for (every_path in site_step:c(length(time_step_1)-1)){
      check <- length(all_shortest_paths(a, every_path, c(every_path+1):numn_nodes, mode = "out")$res)
      if (check==0) {
        site <-0
      }else{
        site <- seq(c(every_path+1),check+every_path,1)
      }
      All_river_paths[every_path,site[1]:site[length(site)]] <- 1
    }
    ST_matrix[temp_connections[site_step],
              c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
    ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+numn_nodes] <- 1
  }else{
    ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+numn_nodes] <- 0
    }
  }
#Lost
if(temp_change==1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes] <- 0
}
#Gain
if(temp_change==-1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+numn_nodes]  <- 0
}

#Temporal indirect links _______________________
#Lost
if(temp_change==1){
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
  for (every_path in site_step:c(length(time_step_1)-1)){
    check <- length(all_shortest_paths(a, every_path+1, c(every_path+1):numn_nodes, mode = "out")$res)
    if (check==0) {
      site <-0
    }else{
      site <- seq(c(every_path+1),check+every_path,1)
    }
    All_river_paths[every_path,site[1]:site[length(site)]] <- 1
  }
  ST_matrix[temp_connections[site_step],
            c(temp_connections[1]+numn_nodes):c(temp_connections[numn_nodes]+numn_nodes)] <- All_river_paths[site_step,]
  
    }
  }
}
ST_matrix_rivers[[river]] <- ST_matrix
}

# Spatiotemporal matrix
ST_matrix_plots <- list()
ST_matrix_out_out <- list()

for (river in 1:length(HOBOS_sites)) {
numn_nodes <- ncol(HOBOS_sites[[river]])-1

out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
for (w in 1:c(length(HOBOS_sites[[river]]$Day)-1)) {

gen_connections <-seq(1,c(length(colnames(HOBOS_sites[[river]]))-1),1)+w*numn_nodes

out <- ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                 gen_connections[1]:gen_connections[length(gen_connections)]]+
       ST_matrix_rivers[[river]][gen_connections[1]:gen_connections[length(gen_connections)],
                 c(gen_connections[1]+numn_nodes):c(gen_connections[length(gen_connections)]+numn_nodes)]
out_out <- out_out+out
}
ST_matrix_out_out[[river]] <- out_out

n<- network(out_out, directed=T, diag=T)

edge_size <- t(out_out)[lower.tri(out_out)]/max(out_out)
edge_col <- rep(0, length(unlist(n$oel)))
for (edg in 1:c(ncol(HOBOS_sites[[river]])-2)) {
row_val <- out_out[edg,][upper.tri(out_out)[edg,]]/max(out_out)
edge_col[n$oel[[edg]]] <- row_val[length(n$oel[[edg]]):1]
}
n %e% "Con_values" <- edge_col

nod_fill_size <- diag(out_out)/max(diag(out_out))
n %v% "Site_values" <- nod_fill_size
ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                 aes(x = x, y = y, xend = xend, yend = yend))+
                            geom_edges(aes(colour=Con_values, size=Con_values/5), arrow=arrow(angle = 20), curvature = 0.15) +
                            geom_nodes(aes(fill=Site_values, size=Site_values*10), color="black" ,shape=21)+
                            scale_color_gradient2(low = "darkred",high = "blue",midpoint = 0.5)+
                            scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
                            theme_classic()+
                            theme(axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  legend.position = "none",
                                  panel.background=element_blank())
}

grid.arrange(ST_matrix_plots[[1]],ST_matrix_plots[[2]],ST_matrix_plots[[3]],
             ST_matrix_plots[[4]],ST_matrix_plots[[5]],ST_matrix_plots[[6]],
             ST_matrix_plots[[7]])


# SpatioTemporal Connectivity________

library(viridis)
ST_connectivity_value <- list()
ST_connectivity_plot <- list()
for (river in 1:length(HOBOS_sites)) {
numn_nodes <- ncol(HOBOS_sites[[river]])-1

ST_matrix<- ST_matrix_rivers[[river]][1:c(numn_nodes*length(HOBOS_sites[[river]]$Day)),1:c(numn_nodes*length(HOBOS_sites[[river]]$Day))]

spt_conn <- c()
leng_correct <- seq(numn_nodes,1,-1)
for (nodes in 1:numn_nodes) {
  spt_conn[nodes] <- sum(apply(ST_matrix,1,sum)[seq(nodes,numn_nodes*length(HOBOS_sites[[river]]$Day)-numn_nodes,numn_nodes)])/leng_correct[nodes]
}
spt_conn<- spt_conn/c(length(HOBOS_sites[[river]]$Day))

ST_connectivity_value[[river]] <- spt_conn  


n<- network(Simple_river_network[[river]], directed=T, diag=T)
  n %v% "CC_values" <- ST_connectivity_value[[river]] #spt_conn
  n %v% "edges_CC_values" <- ST_connectivity_value[[river]]    #spt_conn[1:(numn_nodes-1)]
  
  ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_list[[river]][,4:3]),
                                          aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20)) +
    geom_nodes(aes(fill=CC_values, size=CC_values), color="black" ,shape=21)+
    #scale_color_gradient2()+
    #scale_fill_gradient2()+
    scale_color_viridis(direction = -1)+
    scale_fill_viridis(direction = -1)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())
}

grid.arrange(ST_connectivity_plot[[1]],ST_connectivity_plot[[2]],ST_connectivity_plot[[3]],
             ST_connectivity_plot[[4]],ST_connectivity_plot[[5]],ST_connectivity_plot[[6]],
             ST_connectivity_plot[[7]])


