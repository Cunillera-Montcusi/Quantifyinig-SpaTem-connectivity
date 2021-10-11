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


plot(Sites$Latitud, Sites$Longitud)

VallHorta <- as.data.frame(Sites)%>%
                     filter(Riera=="VH")
VallHorta%>%ggplot(aes(x=Longitud, y=Latitud))+geom_point()+theme_classic()
VallHorta <- VallHorta[2:11,]



detach("package:shp2graph", unload = TRUE)
library(sna)
library(ggnetwork)
#for (e in 1:length(cordenades_xarxes)) {
#  factors <- rep("No_Sampled",nrow(MAPS_xarxes[[e]]))
#  factors[c(nrow(MAPS_xarxes[[e]])-54):nrow(MAPS_xarxes[[e]])] <- "Sampled"
#  CC_values <- PCA_network_results[[e]]

n<- network(out, directed=T, diag=T)
ggplot(n, layout=as.matrix(VallHorta[,4:3]),
                           aes(x = x, y = y, xend = xend, yend = yend))+
    geom_edges(color = "black", size=1, arrow=arrow(angle = 20)) +
    geom_nodes(fill="red", size=5, color="black" ,shape=21)+
    theme_classic()+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          panel.background=element_blank())


# IMPORTANT AQUÏ

ST_matrix_out <- matrix(nrow = ncol(Temp)-1,ncol = ncol(Temp)-1, data=0)
spa_connections <-seq(1,ncol(Temp),1)
for (site_step in 1:c(ncol(Temp)-2)) {
  #Simple spatial links _______________________
  if(time_step_1[site_step]==1){
    ST_matrix_out[spa_connections[site_step],
              spa_connections[site_step]+1] <- 1
  }else{
    ST_matrix_out[spa_connections[site_step],
              spa_connections[site_step]+1] <- 0
  }
}



Temp <- read.csv("Test_1-60days.csv", header = T, sep = ";")
colnames(Temp) <- c("Day","MECHP1","MECHP2","MECHR3","MECHR4","MECHR5","MECHR6","MECHP7","MECHR8","MECHP9","MECHP10")


ST_matrix <- matrix(nrow = 10*length(Temp$Day)+10,ncol = 10*length(Temp$Day)+10, data=0)
for (days in 1:(length(Temp$Day)-1)) {

spa_connections <-seq(1,length(colnames(Temp)),1)+((days-1)*10)

time_step_1 <- Temp[days,2:11]
time_step_2 <- Temp[days+1,2:11]

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
a <- graph.adjacency(ST_matrix[spa_connections[1]:spa_connections[10],
                               spa_connections[1]:spa_connections[10]],mode="directed",diag = FALSE)
All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
for (every_path in 1:c(length(time_step_1)-1)){
  check <- length(all_shortest_paths(a, every_path, c(every_path+1):10, mode = "out")$res)
  if (check==0) {
    site <-0
  }else{
    site <- seq(c(every_path+1),check+every_path,1)
  }
  All_river_paths[every_path,site[1]:site[length(site)]] <- 1
}
ST_matrix[spa_connections[1]:spa_connections[10],
          spa_connections[1]:spa_connections[10]] <- All_river_paths

for (site_step in 1:c(length(time_step_1)-1)) {
#Temporal direct links _______________________
temp_connections <-seq(1,length(colnames(Temp)),1)+((days)*10)  

temp_change <- time_step_1[site_step]-time_step_2[site_step]
#Stable
if(temp_change==0){
  if(time_step_1[site_step]==1){
    #ST_matrix[temp_connections[site_step],
    #          temp_connections[site_step]+10] <- 1
    
    All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
    for (every_path in site_step:c(length(time_step_1)-1)){
      check <- length(all_shortest_paths(a, every_path, c(every_path+1):10, mode = "out")$res)
      if (check==0) {
        site <-0
      }else{
        site <- seq(c(every_path+1),check+every_path,1)
      }
      All_river_paths[every_path,site[1]:site[length(site)]] <- 1
    }
    ST_matrix[spa_connections[site_step],
              c(spa_connections[1]+10):c(spa_connections[10]+10)] <- All_river_paths[site_step,]
  }else{
    ST_matrix[temp_connections[site_step],
              temp_connections[site_step]+10] <- 0
    }
  }
#Lost
if(temp_change==1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+10] <- 0
}
#Gain
if(temp_change==-1){
  ST_matrix[temp_connections[site_step],
            temp_connections[site_step]+10]  <- 0
}

#Temporal indirect links _______________________
#Lost
if(temp_change==1){
  #ST_matrix[temp_connections[site_step],
  #          temp_connections[site_step]+10+1]<- 0.5
  
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
  for (every_path in site_step:c(length(time_step_1)-1)){
    check <- length(all_shortest_paths(a, every_path+1, c(every_path+1):10, mode = "out")$res)
    if (check==0) {
      site <-0
    }else{
      site <- seq(c(every_path+1),check+every_path,1)
    }
    All_river_paths[every_path,site[1]:site[length(site)]] <- 0.5
  }
  ST_matrix[spa_connections[site_step],
            c(spa_connections[1]+10):c(spa_connections[10]+10)] <- All_river_paths[site_step,]
  
    }
  }
}

ST_matrix<- ST_matrix[1:c(10*length(Temp$Day)),1:c(10*length(Temp$Day))]
apply(ST_matrix,1,sum)
plot(apply(ST_matrix,1,sum))

spt_conn <- c(
sum(apply(ST_matrix,1,sum)[seq(11,10*length(Temp$Day)-10,10)])/9,
sum(apply(ST_matrix,1,sum)[seq(12,10*length(Temp$Day)-10,10)])/8,
sum(apply(ST_matrix,1,sum)[seq(13,10*length(Temp$Day)-10,10)])/7,
sum(apply(ST_matrix,1,sum)[seq(14,10*length(Temp$Day)-10,10)])/6,
sum(apply(ST_matrix,1,sum)[seq(15,10*length(Temp$Day)-10,10)])/5,
sum(apply(ST_matrix,1,sum)[seq(16,10*length(Temp$Day)-10,10)])/4,
sum(apply(ST_matrix,1,sum)[seq(17,10*length(Temp$Day)-10,10)])/3,
sum(apply(ST_matrix,1,sum)[seq(18,10*length(Temp$Day)-10,10)])/2,
sum(apply(ST_matrix,1,sum)[seq(19,10*length(Temp$Day)-10,10)])
)/c(length(Temp$Day)-2)

library(viridis)
n<- network(out, directed=T, diag=T)
#n %v% "family" <- factors # Family is an standard name for the categortical variable that we are creating
n %v% "CC_values" <- spt_conn

ggplot(n, layout=as.matrix(VallHorta[,4:3]),
       aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =CC_values), size=spt_conn,arrow=arrow(angle = 20)) +
  geom_nodes(aes(fill=CC_values), size=c(spt_conn*5,0), color="black" ,shape=21)+
  scale_color_gradient2(low = "blue",high = "red",midpoint = 1.5)+
  scale_fill_gradient2(low = "blue",high = "red",midpoint = 1.5)+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.background=element_blank())




