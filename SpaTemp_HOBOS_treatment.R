#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1. Charging & depurating HOBOS Dataset ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")
source("SpaTemp_function.R")

Sites <- read.csv("Raw_HOBOS_Database/Longlat_Rius.csv", header = T, sep = ";", dec = ",")
colnames(Sites) <- c("Riera", "Codi_HOBO","Latitud","Longitud")

#Correction for matching HOBOS -- collecting the  coordinates of the HOBOS that we have  
length(which(as.matrix(dist(Sites[,3:4]))==0))
length(diag(as.matrix(dist(Sites[,3:4]))))
plot(Sites$Longitud, Sites$Latitud)

# We upload the HOBO dataset (1 and 0 defining wet/dry moments)
# We upload the order of the rivers from Upstream to Downstream to arrange the HOBOS into the desired order
## IMPORTANT: This order will also define the direction of the river! 
## The first HOBO in the row of the "Sites_list" must be the first HOBO in the Column of "HOBOS_sites"

HOBOS_sites <- list(
  read.csv("Raw_HOBOS_Database/CA_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/M_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/R_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/SA_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/SC_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/T_HOBOS_data.csv", header = T, sep = ";"),
  read.csv("Raw_HOBOS_Database/VH_HOBOS_data.csv", header = T, sep = ";"))

#Stream order to check
Stream_order <- read.csv("Raw_HOBOS_Database/Stream Order.csv", header = T, sep = ";")

library(dplyr)
# PLotting HOBOS altogether (make the plot window bigger)
par(mfrow=c(3,3))
local_hobos_list <- list()
Sites_list <- list()
for (sit in 1:length(HOBOS_sites)) {
  colnames(HOBOS_sites[[sit]])[1] <- c("Day")
  names_hobos <- colnames(HOBOS_sites[[sit]])[2:length(colnames(HOBOS_sites[[sit]]))]
  local_hobos <- c()
  Sites_list[[sit]] <- Sites%>%filter(Codi_HOBO%in%names_hobos)%>%
                               left_join(Stream_order, by="Codi_HOBO")%>%
                               arrange(UtoD)%>%
                               dplyr::select(Riera,Codi_HOBO,Longitud,Latitud)
#plot(Sites_list[[sit]]$Longitud, Sites_list[[sit]]$Latitud)
}
par(mfrow=c(1,1))

Sites_list

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 2. Dates selection ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

### Autumn dataset - November 2019
BDD <- read.csv2("BiolData/Matriz_otoño.csv", sep=";")
date_correct_Autumn <- min(unique(BDD$Date))
# Edit the dates according the desired format to be entered in the following steps
date_correct_Autumn <- strsplit(date_correct_Autumn,"/") # Split the numbers by the "/"
date_correct_Autumn <- paste(date_correct_Autumn[[1]][3],"-", # Paste the numbers in the correct order and separated by "-" 
                             date_correct_Autumn[[1]][2],"-",
                             date_correct_Autumn[[1]][1],sep = "")

date_correct_Autumn
## All hobos begun at the same date
date_HOBOS <- "2018-07-26"

## All hobos finish at the same date
date_HOBOS_end <- "2019-12-20"

#Beginning
bd <- as.Date(date_HOBOS)

#End
# We change the end for the total dataset (513 days) or for comparing only with Autumn
ed <- as.Date(date_correct_Autumn)

# Difference in number of days will correspond to the number of rows to be selected. 
time_window_beg <- as.numeric(difftime(bd+1, date_HOBOS, units = "days")) 
time_window_beg <- round(time_window_beg)
time_window_end <- as.numeric(difftime(ed, date_HOBOS, units = "days"))
time_window_end <- round(time_window_end)

for (site in 1:length(HOBOS_sites)) {
  HOBOS_sites[[site]] <- HOBOS_sites[[site]][time_window_beg:time_window_end,]
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 3. Distance matrix creation ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Must be a list of distances matrix for each river. 
# Must be a symmetric matrix (upper and lower triangles must be equal). The function already selects which distances 
#are selected in each case.  
# Distances must be small! In this case they are in km
eucl_dist_matrices <- list()
for (river in 1:length(Sites_list)) {
  eucl_dist_matrices[[river]] <-as.matrix(dist(Sites_list[[river]][,4:3],diag = T,upper = T))/1000 
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 4. Network structure creation ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Must be a list with the network strucutre for each river.
# For directed scenarios
Dir_netw_struct_matrices <- list()
for (river in 1:length(HOBOS_sites)) {
  Dir_SpatConnect <- seq(1:c(ncol(HOBOS_sites[[river]])-1))
  netw_struct_matrix <- matrix(nrow =length(Dir_SpatConnect),ncol = length(Dir_SpatConnect),0)
  for (site_step in 1:length(Dir_SpatConnect)-1) {
    netw_struct_matrix[Dir_SpatConnect[site_step],Dir_SpatConnect[site_step]+1] <- 1}
Dir_netw_struct_matrices[[river]] <- netw_struct_matrix
}

# For undirected scenarios 
UNdir_netw_struct_matrices <- list()
for (river in 1:length(HOBOS_sites)) {
  Dir_SpatConnect <- seq(1:c(ncol(HOBOS_sites[[river]])-1))
  netw_struct_matrix <- matrix(nrow =length(Dir_SpatConnect),ncol = length(Dir_SpatConnect),0)
  for (site_step in 1:length(Dir_SpatConnect)) {
    netw_struct_matrix[Dir_SpatConnect[site_step],
                       c(Dir_SpatConnect[1]:length(Dir_SpatConnect))[-site_step]]<- 1}
UNdir_netw_struct_matrices[[river]] <- netw_struct_matrix
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 5. Calculation of ST indices ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Building a directed non weighted network
# Value of link=1
# Value of NO link=0
Dir_NonW_Net <- spat_temp_index(Inermitence_dataset = HOBOS_sites,
                                Sites_coordinates=Sites_list,
                                direction="directed", sense = "out",
                                Network_stru =Dir_netw_struct_matrices,
                                weighting=FALSE,dist_matrices = NULL,
                                weighting_links =FALSE,link_weights = NULL,
                                legacy_effect =1, legacy_lenght = 1,
                                value_S_LINK=1,
                                value_T_LINK=1,
                                value_NO_S_link=0,
                                value_NO_T_link=0,
                                Network_variables=F,print.plots=F,print.directory="Figure/")

# Building a directed Weighted network
# Value of link= 0.1
# Value of NO link=1
# - I am evaluationg "resistance" to move from one to another. 
# SO, higher numbers mean High resistance
Dir_WEIG_Net <- spat_temp_index(Inermitence_dataset =HOBOS_sites,
                                Sites_coordinates=Sites_list,
                                direction="directed", sense = "out",
                                weighting=TRUE, dist_matrices = eucl_dist_matrices,
                                Network_stru =Dir_netw_struct_matrices,
                                weighting_links =FALSE,link_weights = NULL,
                                legacy_effect =1, legacy_lenght = 1,
                                value_S_LINK=0.1,
                                value_T_LINK=0.1,
                                value_NO_S_link=1,
                                value_NO_T_link=1,Virid_option = "B",
                                Network_variables=F,print.plots=F,print.directory="Figure/")

# Building a Undirected non weighted network
# Value of link=1
# Value of NO link=0
UnD_NonW_Net <- spat_temp_index(Inermitence_dataset =HOBOS_sites,
                                Sites_coordinates=Sites_list,
                                direction="undirected", sense = "all",
                                weighting=FALSE,dist_matrices = NULL,
                                Network_stru =UNdir_netw_struct_matrices, 
                                weighting_links =FALSE,link_weights = NULL,
                                legacy_effect =1, legacy_lenght = 1,
                                value_S_LINK=1,
                                value_T_LINK=1,
                                value_NO_S_link=0,
                                value_NO_T_link=0,
                                Network_variables=F,print.plots=F,print.directory="Figure/")

# Building a Undirected Weighted network
# Value of link= 0.1
# Value of NO link=1
# - I am evaluationg "resistance" to move from one to another. 
# SO, higher numbers mean High resistance
UnD_WEIG_Net <- spat_temp_index(Inermitence_dataset = HOBOS_sites,
                                Sites_coordinates=Sites_list,
                                direction="undirected",  sense = "all",
                                weighting=TRUE, dist_matrices = eucl_dist_matrices,
                                Network_stru =UNdir_netw_struct_matrices,
                                weighting_links =FALSE,link_weights = NULL,
                                legacy_effect =1, legacy_lenght = 1,
                                value_S_LINK=0.1,
                                value_T_LINK=0.1,
                                value_NO_S_link=1,
                                value_NO_T_link=1,Virid_option = "B",
                                Network_variables=F,print.plots=F,print.directory="Figure/")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 5. Data extraction and indices treatment ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# 5.1. Directed NONWEIGHTED river network OUTPUTS ####
# STcon 
NonW_ST_connectivity_value <- Dir_NonW_Net$STcon
NonW_ST_directed_output <- data.frame(NonW_Dir_con=unlist(NonW_ST_connectivity_value))
# STconmat 
NonW_ST_matrix_out_out <- Dir_NonW_Net$STconmat


# 5.2. Directed WEIGHTED river network OUTPUTS ####
# STcon
WEIG_ST_connectivity_value <- Dir_WEIG_Net$STcon
WEIG_ST_directed_output <- data.frame(WEIG_Dir_con=unlist(WEIG_ST_connectivity_value))
# STconmat 
WEIG_ST_matrix_out_out <- Dir_WEIG_Net$STconmat

# 5.3. Undirected river network OUTPUTS ####
# STcon 
Un_NonW_ST_connectivity_value <- UnD_NonW_Net$STcon
Un_NonW_ST_output <- data.frame(Un_NonW_con=unlist(Un_NonW_ST_connectivity_value))
# STconmat
Un_NonW_ST_matrix_out_out <- UnD_NonW_Net$STconmat

# 5.4 Undirected WEIGHTED river network OUTPUTS ####
# STcon
Un_WEIG_ST_connectivity_value <- UnD_WEIG_Net$STcon
Un_WEIG_ST_output <- data.frame(Un_WEIG_con=unlist(Un_WEIG_ST_connectivity_value))
# STconmat 
Un_WEIG_ST_matrix_out_out <- UnD_WEIG_Net$STconmat


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 6. Final values__________________________________________________________________ ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# We add the stream ID and the order (UtoD) to each output to later analyse them. 
# We retrieve this information from the "Stream_order" dataset where each HOBO code, Stram and its position on the 
#upstream downstream direction were registered
# This two following vectors are created for being used in other scripts
HOB_riv_ID<- Stream_order$ï..Riera
ups_dos <- Stream_order$UtoD

NonW_ST_directed_out <- data.frame(ID=Stream_order$ï..Riera,DtoU=Stream_order$UtoD,
                                   NonW_ST_directed_output)%>%
                        mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))

WEIG_ST_directed_out <- data.frame(ID=Stream_order$ï..Riera,DtoU=Stream_order$UtoD,
                                   WEIG_ST_directed_output)%>%
                        mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))

Un_NonW_ST_out <- data.frame(ID=Stream_order$ï..Riera,DtoU=Stream_order$UtoD,
                             Un_NonW_ST_output)%>%
                        mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))

Un_WEIG_ST_out <- data.frame(ID=Stream_order$ï..Riera,DtoU=Stream_order$UtoD,
                             Un_WEIG_ST_output)%>%
                        mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))

# Final values HOBOS dataframe
NonW_ST_matrix_out_out
WEIG_ST_matrix_out_out
Un_NonW_ST_matrix_out_out
Un_WEIG_ST_matrix_out_out

NonW_ST_directed_out
WEIG_ST_directed_out
Un_NonW_ST_out
Un_WEIG_ST_out

# Result printing outputs
names_scenarios <- c("Scenario A1",
                     "Scenario A2",
                     "Scenario B1",
                     "Scenario B2")
for (n in 1:length(NonW_ST_matrix_out_out)){
#Scenario A1 STconmat 
write.table(rbind(
            cbind(Scenario=names_scenarios[1],Sites_list[[n]], round(NonW_ST_matrix_out_out[[n]],3)),
            cbind(Scenario=names_scenarios[2],Sites_list[[n]], round(WEIG_ST_matrix_out_out[[n]],3)),
            cbind(Scenario=names_scenarios[3],Sites_list[[n]], round(Un_NonW_ST_matrix_out_out[[n]],3)),
            cbind(Scenario=names_scenarios[4],Sites_list[[n]], round(Un_WEIG_ST_matrix_out_out[[n]],3))),
            paste("Table_Results/",unique(Sites_list[[n]][,1]),"stream", "STconmat.txt",sep = "_"), sep = ",")
}

#Scenario A1 STcon 
write.table(cbind(rbind(Sites_list[[1]],Sites_list[[2]],Sites_list[[3]],
                  Sites_list[[4]],Sites_list[[5]],Sites_list[[6]],Sites_list[[7]]),
                  NonW_ST_directed_out[,1:2],
                  round(NonW_ST_directed_out[,3:6],3)), "Table_Results/All streams_STcon_A1.txt", sep = ",")
#Scenario A2 STcon 
write.table(cbind(rbind(Sites_list[[1]],Sites_list[[2]],Sites_list[[3]],
                        Sites_list[[4]],Sites_list[[5]],Sites_list[[6]],Sites_list[[7]]),
                  WEIG_ST_directed_out[,1:2],
                  round(WEIG_ST_directed_out[,3:6],3)), "Table_Results/All streams_STcon_A2.txt", sep = ",")
#Scenario B1 STcon 
write.table(cbind(rbind(Sites_list[[1]],Sites_list[[2]],Sites_list[[3]],
                        Sites_list[[4]],Sites_list[[5]],Sites_list[[6]],Sites_list[[7]]),
                  Un_NonW_ST_out[,1:2],
                  round(Un_NonW_ST_out[,3:6],3)), "Table_Results/All streams_STcon_B1.txt", sep = ",")
#Scenario B2 STcon 
write.table(cbind(rbind(Sites_list[[1]],Sites_list[[2]],Sites_list[[3]],
                        Sites_list[[4]],Sites_list[[5]],Sites_list[[6]],Sites_list[[7]]),
                  Un_WEIG_ST_out[,1:2],
                  round(Un_WEIG_ST_out[,3:6],3)), "Table_Results/All streams_STcon_B2.txt", sep = ",")


