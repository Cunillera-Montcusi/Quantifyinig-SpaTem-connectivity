#__________________________________________________________________________________________________________#
#__________________________________________________________________________________________________________#
#______________ spat_temp_index function __________________________________________________________________#
#__________________________________________________________________________________________________________#
#__________________________________________________________________________________________________________#

# We present the function "spat_temp_index" a function to calculate Spatiotemporal connectiviy indices based on 
#spatiotemporal graph networks. The aim of this function is to provide a methodological framework from which to 
#calculate these indices based on high-frequency (or frequency-based) information obtained from natural ecosystems. 

# The current script as well as all the information contained in it and its tutorial are part of the article titled: 
#Navigating through space and time: a methodological approach to quantify spatiotemporal connectivity using flow intermittence data as a study case.
#David Cunillera-Montcusi1,2,3,4*, Jose Maria Fernandez-Calero1,2, Sebastian Polsterl5, Julia Valera1, Roger Argelich1, Nuria Cid1,6, Nuria Bonada1,2,
#Miguel Canedo-Arguelles1,7,8
#1- FEHM-Lab (Freshwater Ecology, Hydrology and Management), Departament de Biologia Evolutiva, Ecologia i Ciencies Ambientals, Facultat de Biologia,Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Spain.
#2- Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Spain.
#3. Departamento de Ecologia y Gestion Ambiental, Centro Universitario Regional del Este (CURE), Universidad de la Repulica, Tacuarembo? s/n, Maldonado, Uruguay.
#4. GRECO, Institute of Aquatic Ecology, University of Girona, Girona, Spain
#5- The Lab for Artificial Intelligence in Medical Imaging (AI-Med), Department of Child and Adolescent Psychiatry, Ludwig-Maximilians-Universitat, Waltherstrae 23, 80337 Munich, Germany
#6- IRTA  Marine and Continental Waters Programme, Ctra de Poble Nou Km 5.5, E43540, La Rapita, Catalonia, Spain 
#7- Institut de Recerca de l'Aigua (IdRA), Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Catalonia, Spain.
#8- Institute of Environmental Assessment and Water Research (IDAEA-CSIC), Carrer de Jordi Girona, 18-26, 08034 Barcelona

#*Corresponding author: david.cunillera@dcm.cat

# The script and the follwouing functions have been written by David Cunillera-Montcusi. For any question, comment or 
#feedback contact him at: david.cunillera@dcm.cat 

# Mainly, the information that should be contained in the datasets must represent some sort of habitat availability 
# like habitat presence/absence (e.g. aquatic habitats wet/dry phases) or habitat connectivity strengths (e.g. flows 
# in rivers) or any other type of interpretation that has ecological meaning and that can be related with Spatial and 
# temporal patterns. 

# ONE main information must be provided and must be considered before to use the function: THE SPATIOTEMPORAL INFORMATION AND STRUCTURE
#The function is thought for streams/rivers which have a stron directional pattern (Upstream to downstream). Therefore, the table provided
#must respect this directionality as well as the required format to let the method work properly. 

# See below, an example of the type of matrix that must be entered in the function as HOBOS_dataset

#data.frame(
#MonitoredDays= c("Day1","Day2","Day3","Day4"), # An identifier for the monitored days from ORDERED from the "oldest" to the "newest"
#StreamSite1=c(0,1,1,1), # More upstream site water presence record (1= water presence, 0= water absence)
#StreamSite2=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite3=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite4=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite5=c(0,1,0,1), # Following monitored site water presence record in descending order 
#StreamSite6=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite7=c(0,1,0,1), # Following monitored site water presence record in descending order 
#StreamSite8=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite9=c(1,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite10=c(1,1,1,1),# Following monitored site water presence record in descending order 
#)

# This information is KEY to preserve the structure and functioning of the function. The system does not necessarily need to quantify 
#water absence/presence but each cell value must represent a feature defining connectivity "on" or "off" and that can be transmitted 
#to built ecologically meaningful links in a spatiotemporal graph. 
# The meaning of a "yes" (Input matrix value= 1) or a "no" (Input matrix value= 0) can be later modified by the values assigned to  
# the "LINK" (spatial and temporal) within the figure. 

# direction either "directed" or "undirected" --> Undirected considers that all neighbours are equally reachable. 

# weighting either FALSE or TRUE --> The value of the weight is defined by the matrix added (must be in a distance matrix format) and 
#can consider any type of distance between pairs of monitored sites (euclidean, environmental, topographic, ...)
#dist_matrices--> corresponds to the attached matrix representing the "distances" between sites. 

# value_LINK is the value attributed for each "effective link", which is a connection between two nodes (a line)
# - value_S_LINK for spatial links
# - value_T_LINK for temporal links
# value_NO_link is the value for each "effective disconnection", which is a "connection"void" between (a white space)
# - value_NO_S_LINK for spatial links
# - value_NO_T_LINK for temporal links

# WARNING ! All values must be entered in a list format, even if only 1 matrix is used. This is done in order to allow the possibility 
# to enter several streams in one call and obtain the results according to that.   

spat_temp_index <- function(HOBOS_dataset, 
                            Sites_coordinates,
                            direction,
                            sense,
                            weighting,
                            Network_stru,
                            dist_matrices,
                            value_S_LINK=1,
                            value_T_LINK=1,
                            value_NO_S_link=0,
                            value_NO_T_link=0,
                            Network_variables=TRUE,
                            print.plots=TRUE,
                            print.directory){
  
  require(gridExtra);require(igraph);require(shp2graph);require(sna,quietly = T,warn.conflicts = F);
  require(tidyverse);require(viridis);require(doParallel);require(ggnetwork)
  
  if(print.plots==TRUE){
    if(is.character(print.directory)==T){
      cat("Hope that the plot.directory is the directory where you want to save the photos with a last '/' on it ","\n")}
    else{
      return("You must write print.directory as character and make sure that the directory is the correct folder","\n")}
  }else{
    cat("Plots will not be printed in a folder but created together with the outputs","\n")
  }
  
  if(direction=="directed"){
    cat("Your river will be considered as a directed graph","\n")}
  if(direction=="undirected"){
    cat("Your river will be considered as an undirected graph","\n")}
  
  # We select the corresponding distance matrix
  if(weighting==T){cat("Your connectivity will be Weighted, connections will be multiplied by 'dist_matrices'","\n")}
  if(weighting==F){cat("Your connectivity will be NON weighted, connections will not be multiplied by any distance matrix","\n")}
  
  # Simple river network________________####
  # Building river networks based on just directional network.
  Simple_river_network <- list()
  Simple_river_network_maps <- list()
  for (river in 1:length(HOBOS_dataset)) {
    ST_matrix_out <- matrix(nrow = ncol(HOBOS_dataset[[river]])-1,ncol = ncol(HOBOS_dataset[[river]])-1, data=0)
    spa_connections <-seq(1,ncol(HOBOS_dataset[[river]])-1,1)
    time_step_1 <- rep(1,ncol(HOBOS_dataset[[river]])-1)
    
    for (site_step in 1:c(ncol(HOBOS_dataset[[river]])-1)) {
      #Simple spatial links _______________________
      if(time_step_1[site_step]==1){
        ST_matrix_out[spa_connections[site_step],which(Network_stru[[river]][site_step,]==1)] <- 1
      }else{
        ST_matrix_out[spa_connections[site_step],which(Network_stru[[river]][site_step,]==1)] <- 0
      }
    }
    Simple_river_network[[river]] <- ST_matrix_out
    
    require(sna,quietly = T,warn.conflicts = F)
    n<- network(Simple_river_network[[river]], directed=T, diag=T)
    Simple_river_network_maps[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
                                                 aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(color = "black", size=1, arrow=arrow(angle = 20)) +
      geom_nodes(fill="red", size=5, color="black" ,shape=21)+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background=element_blank())
  }
  ####_______________________________________________________________________
  if(print.plots==TRUE){
    png(filename = paste(print.directory,"Simple_river_network_maps.png"),
        width = 715*6,height = 448*6,units = "px",res = 300)
    grid.arrange(arrangeGrob(grobs = Simple_river_network_maps,ncol=ceiling(length(HOBOS_dataset)/2)),
                 top="Simple river network")
    dev.off()}
  
  ####_______________________________________________________________________
  # River network ####
  ####_______________________________________________________________________
  # River network matrix BUILDING ####
  ####_______________________________________________________________________
  pack_check <- search()
  pack_check_val <- length(which(pack_check=="package:sna"))
  if(pack_check_val>0){detach("package:sna", unload = TRUE)}
  # Below there is the function who builds the MATRIX ponderating SPATIAL lINKS=1 and TEMPORAL LINKS=1
  ST_matrix_rivers <- list()
  ST_directed_Ocloseness_rivers <- list()
  ST_directed_Allcloseness_rivers <- list()
  ST_directed_betweennes_rivers <- list()
  
  #Parallelization parameters
  registerDoParallel(cores = detectCores()-1)
  
  out_Matrix_LIST <- list()
  pack_check <- search()
  pack_check_val <- length(which(pack_check=="package:sna"))
  if(pack_check_val>0){detach("package:sna", unload = TRUE)}
  out_Matrix_LIST <- foreach(river=1:length(HOBOS_dataset))%dopar%{
    #for (river in 1:length(HOBOS_dataset)) { - With this it takes 6'26''
    # We calculate the number of nodes of our network (used along the function)  
    numn_nodes <- ncol(HOBOS_dataset[[river]])-1
    
    if(weighting==TRUE){dist_matr <- dist_matrices[[river]]}
    
    # We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
    ### This matrix is the "giant" themplate where we will put all the values.
    ST_matrix <- matrix(nrow = numn_nodes,ncol = numn_nodes*2, data=0)
    ST_matrix_netwGraph <- matrix(nrow = numn_nodes,ncol = numn_nodes, data=0)
    
    ST_Oclosenness_matrix <- matrix(length(HOBOS_dataset[[river]][,1]),
                                    numn_nodes, data=0)
    ST_Allclosenness_matrix <- matrix(length(HOBOS_dataset[[river]][,1]),
                                      numn_nodes, data=0)
    ST_betweennes_matrix <- matrix(length(HOBOS_dataset[[river]][,1]),
                                   numn_nodes, data=0)
    
    # Once created the template we start to fill it for every day
    ### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
    for (days in 1:c(length(HOBOS_dataset[[river]][,1])-1)) {
      print(paste("We are at time unit", days, "of", (length(HOBOS_dataset[[river]][,1])-1)))
      # First we define the spatial connections of the matrix
      ### Also known as the rows or columns at which we have to add the values of the connections 
      spa_connections <-seq(1,length(colnames(HOBOS_dataset[[river]]))-1,1)#+((days-1)*numn_nodes)
      
      # We obtain the time steps:
      ## time_step_1 is the present
      ## time_step_2 is the following step (the close future)
      time_step_1 <- HOBOS_dataset[[river]][days,2:ncol(HOBOS_dataset[[river]])]
      time_step_2 <- HOBOS_dataset[[river]][days+1,2:ncol(HOBOS_dataset[[river]])]
      
      #Simple fluvial network_______________________
      ## This step fills "the diagonal" of each time_step following the direction of the river
      ## it basically connects the river in a dendritic structure.
      for (site_step in 1:c(length(time_step_1))) {
        if(time_step_1[site_step]==1){
          ST_matrix_netwGraph[spa_connections[site_step],
                              c(spa_connections[1]:spa_connections[numn_nodes])] <- as.numeric(Network_stru[[river]][site_step,])
          # We weight
          if(weighting==T){ST_matrix_netwGraph[spa_connections[site_step],
                                               c(spa_connections[1]:spa_connections[numn_nodes])] <-ST_matrix_netwGraph[spa_connections[site_step],
                                                                                                                        c(spa_connections[1]:spa_connections[numn_nodes])]*dist_matr[site_step,]}
        }else{
          ST_matrix_netwGraph[spa_connections[site_step],
                              c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 0
          # We weight
          if(weighting==T){ST_matrix_netwGraph[spa_connections[site_step],
                                               c(spa_connections[1]:spa_connections[numn_nodes])] <-ST_matrix_netwGraph[spa_connections[site_step],
                                                                                                                        c(spa_connections[1]:spa_connections[numn_nodes])]*dist_matr[site_step,]}
        }
      }
      
      # FLuvial SPATIAL links ___________________________________________________________________________________________________________________
      # Now the party begins. 
      ## Here we fill the matrix section corresponding to the time_step based on the river graph based on a dendritic. 
      require(igraph)
      # We create the graph
      a <- graph.adjacency(ST_matrix_netwGraph[spa_connections[1]:spa_connections[numn_nodes],
                                               spa_connections[1]:spa_connections[numn_nodes]], 
                           mode=direction,diag = FALSE)
      
      if(Network_variables==T){
        ST_Oclosenness_matrix[days,] <- closeness(a, mode = sense,normalized = T)
        ST_Oclosenness_matrix[days,which(time_step_1==0)] <- 0
        ST_Allclosenness_matrix[days,] <- closeness(a, mode = "all",normalized = T)
        ST_betweennes_matrix[days,] <- betweenness(a)
      }
      
      # We create the matrix where we will drop the information of the shortest paths.
      ## We will fill "1" or "0" according to the shortest paths. 
      All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_S_link)
      #All_river_paths[upper.tri(All_river_paths)] <- value_NO_S_link
      
      # For each path (e.g., from node 1 to node 7) we "check" the length of the shortest path. 
      ## check = 0 means that the graph is disconnected.
      ## check bigger than 0 means that the graph is connected.
      for (every_path in 1:c(length(time_step_1))){
        check <- length(all_shortest_paths(a, every_path, 1:numn_nodes, mode = sense)$res)
        if (check==0) {
          site <-0
        }else{
          # If bigger than 0. we create a sequence from the path to downstream.
          neigh <- (c(1:numn_nodes)[-every_path])
          connect_loc <- all_shortest_paths(a, every_path,neigh, mode = sense)$nrgeo
          connect_loc[every_path] <- 0
          site <- which(connect_loc>0)
        }
        # We fill the "All_river_paths" with 1 on the connections concerning to each "row" or node.
        ## Site is the vector with the connections (follwing the river downstream).
        ## When "0" site does not correspond to any row... so the "1" does not go anywhere. 
        All_river_paths[every_path,site] <- value_S_LINK
        # We weight
        if(weighting==T){All_river_paths[every_path,] <- All_river_paths[every_path,]*dist_matr[every_path,]}
      }
      
      # We add the "All_river_paths" filled for each node in the "big" matrix specific sites
      ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths+ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                                                                                             spa_connections[1]:spa_connections[numn_nodes]]
      
      # In the following lines we continue the party towards temporal steps
      for (site_step in 1:length(time_step_1)) {
        # We created "All_river_paths" for temporal
        All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_NO_T_link)
        
        # FLuvial TEMPORAL DIRECT links ___________________________________________________________________________________________________________________
        ## We generate the temporal connectins
        temp_connections <-seq(1+numn_nodes,length(colnames(HOBOS_dataset[[river]]))-1+numn_nodes,1)#+((days)*numn_nodes) 
        
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
            # We then do the same as before, check, substitute and add "1" or 0 depending if the connection 
            # following the river is flow.
            #for (every_path in 1:length(time_step_1)){
            check <- length(all_shortest_paths(a, site_step, 1:numn_nodes, mode = sense)$res)
            if (check==0) {
              site <-0
            }else{
              neigh <- (c(1:numn_nodes)[-site_step])
              connect_loc <- all_shortest_paths(a, site_step,neigh, mode = sense)$nrgeo
              connect_loc[site_step] <- 0
              site <- which(connect_loc>0)
            }
            All_river_paths[site_step,site] <- value_T_LINK
            # We weight
            if(weighting==T){All_river_paths[site_step,] <- All_river_paths[site_step,]*dist_matr[site_step,]}
            
            # TEMPORAL LINKS are filled in the "future" of our current matrix. This means that we are filling the matrix in 
            # in the diagonal of our "time step" for spatial links but we add the temporal links in the following time step. 
            # so, we evaluate here the present (time step 1) and the future (time step 2) but we register it as the past of the future (at time step 2)
            ST_matrix[spa_connections[site_step],
                      temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths[site_step,]+ST_matrix[spa_connections[site_step],
                                                                                                                 temp_connections[1]:temp_connections[numn_nodes]]
            # Here we add the temporal "link" between "himself". If the link is stable and connected (from 1 to 1), we fill the 
            # diagonal value accordingly. Therefore, we will be able to evaluate the relationship between "himself". Kind of 
            # Tot_Num indicator.  
            ST_matrix[spa_connections[site_step],
                      temp_connections[site_step]] <- value_T_LINK+ST_matrix[spa_connections[site_step],temp_connections[site_step]]
          }else{# Here we check if the temporal change implies going from 0 to 0 (so a stable disconnected link). Then we put 0
            # We weight
            if(weighting==T){All_river_paths[site_step,] <- All_river_paths[site_step,]*dist_matr[site_step,]}
            ST_matrix[spa_connections[site_step],
                      temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths[site_step,]+ST_matrix[spa_connections[site_step],
                                                                                                                 temp_connections[1]:temp_connections[numn_nodes]]
          }
        }
        
        #Lost links (when temp_change=1)
        ## This just needs to be filled with zeros... so no need to use "All_river_paths"
        if(temp_change==1){
          if(weighting==T){All_river_paths[site_step,] <- All_river_paths[site_step,]*dist_matr[site_step,]}
          ST_matrix[spa_connections[site_step],
                    temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths[site_step,]+ST_matrix[spa_connections[site_step],
                                                                                                               temp_connections[1]:temp_connections[numn_nodes]]
        }
        #Gained links (when temp_change=-1)
        ## It is a "gain" but it means that "in the present" (time step 1), the node is still disconnected. So it =0
        if(temp_change==-1){
          if(weighting==T){All_river_paths[site_step,]*dist_matr[site_step,]}
          ST_matrix[spa_connections[site_step],
                    temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths[site_step,]+ST_matrix[spa_connections[site_step],
                                                                                                               temp_connections[1]:temp_connections[numn_nodes]]
        }
        
        # FLuvial TEMPORAL INDIRECT links ___________________________________________________________________________________________________________________
        ## Indirect links are those links defined here as the ones that are "lost for the first time" (so temp_change=1).
        ## When this occurs we asign an "extra" 1 in that particular case. Assuming that when the node dries for the first time
        ## there is an increase in "dispersal" (downstream directed).
        ### This only occurs when there is a loss of a previously wet node (from 1 in the present to 0 in the future).
        if(temp_change==1){
          #All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
          #All_river_paths[upper.tri(All_river_paths)] <- value_NO_T_link
          check <- length(all_shortest_paths(a, site_step, 1:numn_nodes, mode = sense)$res)
          if (check==0) {
            site <-0
          }else{
            neigh <- (c(1:numn_nodes)[-site_step])
            connect_loc <- all_shortest_paths(a, site_step,neigh, mode = sense)$nrgeo
            connect_loc[site_step] <- 0
            site <- which(connect_loc>0)
          }
          # We fill the sites with the value
          All_river_paths[site_step,site] <- value_T_LINK
          # We weight
          if(weighting==T){All_river_paths[site_step,] <- All_river_paths[site_step,]*dist_matr[site_step,]}
          # We pass it to the main matrix
          ST_matrix[spa_connections[site_step],
                    temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths[site_step,]+ST_matrix[spa_connections[site_step],
                                                                                                               temp_connections[1]:temp_connections[numn_nodes]]
        }
        
      }# Site_step closing
    }# Days closing
    
    out_Matrix <- list(ST_matrix,ST_Oclosenness_matrix,ST_Allclosenness_matrix,ST_betweennes_matrix)
    out_Matrix_LIST[[river]] <- out_Matrix
  }# Loop for every river entered in the lists
  
  # Exctracring the results into different lists
  for (river in 1:length(HOBOS_dataset)) {ST_matrix_rivers[[river]] <- out_Matrix_LIST[[river]][[1]]}
  
  if(Network_variables==T){
    for (river in 1:length(HOBOS_dataset)) {
      ST_directed_Ocloseness_rivers[[river]] <- out_Matrix_LIST[[river]][[2]]
      ST_directed_Allcloseness_rivers[[river]] <- out_Matrix_LIST[[river]][[3]]
      ST_directed_betweennes_rivers[[river]] <- out_Matrix_LIST[[river]][[4]]
    }
  }
  
  ####_______________________________________________________________________
  # STconmat calculaiton ####
  ####_______________________________________________________________________
  # Find below the lines to calculate the "collapsing" matrix that just summs all the values of all the SPATIOTEMPORAL matrix
  # These pairwise matrix is called the STconmat.   
  ST_matrix_plots <- list()
  ST_matrix_out_out <- list()
  
  for (river in 1:length(HOBOS_dataset)) {
    numn_nodes <- ncol(HOBOS_dataset[[river]])-1
    temp_connections <-seq(1+numn_nodes,length(colnames(HOBOS_dataset[[river]]))-1+numn_nodes,1)
    spa_connections <-seq(1,length(colnames(HOBOS_dataset[[river]]))-1,1)
    
    # We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
    out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
    
    Spatial_matrix <-  ST_matrix_rivers[[river]][,spa_connections]
    Temporal_matrix <- ST_matrix_rivers[[river]][,temp_connections]
    
    out_out <- Spatial_matrix+Temporal_matrix
    out_out <- out_out/c(length(HOBOS_dataset[[river]][,1])-1)
    
    # We save the collapsed matrix
    ST_matrix_out_out[[river]] <- out_out
    
    # Following lines are just a plotting schema to obtain the graphic representation. 
    # Note that the values are "scaled" to one for colours and sizes.
    n<- network(out_out, directed=T, diag=F)
    
    edge_col <- rep(0, length(unlist(n$oel)))
    
    for (edg in 1:c(ncol(HOBOS_dataset[[river]])-2)) {
      no_diag_out_out <- out_out
      diag(no_diag_out_out) <- 0
      # This line is a bit tricky but key for coloring properly the edges according to their value! 
      edge_col[n$oel[[edg]]] <- rev(c(no_diag_out_out[edg,]/max(no_diag_out_out))[-seq(edg,1)])
    }
    
    n %e% "Con_values" <- edge_col
    n %e% "Con_values_SIZE" <- ((edge_col)^5)*2
    nod_fill_size <- diag(out_out)/max(diag(out_out))
    n %v% "Site_values" <- nod_fill_size
    ST_matrix_plots[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
                                       aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(aes(colour=Con_values, size=Con_values_SIZE), alpha=0.2,
                 arrow=arrow(angle = 20), curvature = 0.15) +
      #geom_nodes(aes(fill=Site_values, size=Site_values*10), color="black" ,shape=21)+
      geom_nodes(size=3, fill="grey40", color="black" ,shape=21)+
      scale_color_viridis(direction = -1)+
      scale_fill_viridis(direction = -1)+
      #scale_color_gradient2(low = "brown",high = "blue",midpoint = 0.5)+
      #scale_fill_gradient2(low ="darkred",high = "darkgreen",midpoint = 0.5)+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background=element_blank())
  }
  ####_______________________________________________________________________
  if(print.plots==TRUE){
    png(filename = paste(print.directory,"STconmat.png"),
        width = 715*6,height = 448*6,units = "px",res = 300)
    grid.arrange(arrangeGrob(grobs = ST_matrix_plots,ncol=ceiling(length(HOBOS_dataset)/2)),
                 top="ST Dir NonW connectivity Matrix")
    dev.off()}
  
  ####_______________________________________________________________________
  # STcon calculaiton ####
  ####_______________________________________________________________________
  ST_connectivity_value <- list()
  ST_connectivity_plot <- list()
  for (river in 1:length(HOBOS_dataset)) {
    # We already know this value
    numn_nodes <- ncol(HOBOS_dataset[[river]])-1
    temp_connections <-seq(1+numn_nodes,length(colnames(HOBOS_dataset[[river]]))-1+numn_nodes,1)
    spa_connections <-seq(1,length(colnames(HOBOS_dataset[[river]]))-1,1)
    
    # We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
    out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)
    
    Spatial_matrix <-  ST_matrix_rivers[[river]][,spa_connections]
    Temporal_matrix <- ST_matrix_rivers[[river]][,temp_connections]
    
    out_out <- Spatial_matrix+Temporal_matrix
    
    # "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
    # considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
    leng_correct <- c()
    aa <- graph_from_adjacency_matrix(as.matrix(Network_stru[[river]]),mode = "directed")
    for (neigg in 1:numn_nodes) {
      neighbour<- all_shortest_paths(aa, from = neigg, to = c(1:numn_nodes)[-neigg],mode = sense)$nrgeo
      neighbour[neigg] <- 0
      leng_correct[neigg] <- length(which(neighbour>0))
    }
    
    spt_conn <-apply(out_out,1,sum)/leng_correct
    # We divide by the number of days so we obtain the "per day" values
    spt_conn<- spt_conn/c(length(HOBOS_dataset[[river]][,1])-1)
    
    ST_connectivity_value[[river]] <- spt_conn  
    
    n<- network(as.matrix(Network_stru[[river]]), directed=T, diag=T)
    n %v% "CC_values" <- c(ST_connectivity_value[[river]][1:numn_nodes]) 
    n %v% "edges_CC_values" <- c(ST_connectivity_value[[river]][1:numn_nodes]) 
    
    ST_connectivity_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
                                            aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(aes(color =edges_CC_values, size=edges_CC_values),arrow=arrow(angle = 20),curvature = 0.15) +
      #geom_nodes(aes(size=CC_values), color="black" ,shape=16)+
      geom_nodes(size=3, fill="grey40", color="black" ,shape=21)+
      scale_color_viridis(direction = -1)+
      scale_fill_viridis(direction = -1)+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.background=element_blank())
  }
  
  if(print.plots==TRUE){
    png(filename = paste(print.directory,"STcon.png"),
        width = 715*6,height = 448*6,units = "px",res = 400)
    grid.arrange(arrangeGrob(grobs = ST_connectivity_plot,ncol=ceiling(length(HOBOS_dataset)/2)),
                 top="ST Dir NonW connectivity")
    dev.off()}
  
  if(Network_variables==T){
    ####_______________________________________________________________________
    # ST Out closeness ####
    ####_______________________________________________________________________
    ST_Oclo_plot <- list()
    for (river in 1:length(HOBOS_dataset)) {
      mean_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,mean)
      sd_Oclo <- apply(ST_directed_Ocloseness_rivers[[river]],2,sd)
      
      n<- network(Simple_river_network[[river]], directed=T, diag=T)
      n %v% "CC_values" <- mean_Oclo
      n %v% "Sd_values" <- mean_Oclo+sd_Oclo
      n %v% "edges_CC_values" <- mean_Oclo[1:length(mean_Oclo)-1]
      
      ST_Oclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
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
    
    if(print.plots==TRUE){
      png(filename = paste(print.directory,"SToutclosennes.png"),
          width = 715*6,height = 448*6,units = "px",res = 400)
      grid.arrange(arrangeGrob(grobs = ST_Oclo_plot,ncol=ceiling(length(HOBOS_dataset)/2)),
                   top="ST Dir NonW Out closennes")
      dev.off()}
    
    
    ####_______________________________________________________________________
    # ST closeness ####
    ####_______________________________________________________________________
    ST_Allclo_plot <- list()
    for (river in 1:length(HOBOS_dataset)) {
      mean_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,mean)
      sd_Allclo <- apply(ST_directed_Allcloseness_rivers[[river]],2,sd)
      
      n<- network(Simple_river_network[[river]], directed=T, diag=T)
      n %v% "CC_values" <- mean_Allclo
      n %v% "Sd_values" <- mean_Allclo+sd_Allclo
      n %v% "edges_CC_values" <- mean_Allclo[1:length(mean_Allclo)-1]
      
      ST_Allclo_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
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
    
    if(print.plots==TRUE){
      png(filename = paste(print.directory,"STclosennes.png"),
          width = 715*6,height = 448*6,units = "px",res = 400)
      grid.arrange(arrangeGrob(grobs = ST_Allclo_plot,ncol=ceiling(length(HOBOS_dataset)/2)),
                   top="ST Dir NonW All closennes")
      dev.off()}
    
    
    ####_______________________________________________________________________
    # ST Betweenness calculaiton ####
    ####_______________________________________________________________________
    ST_betw_plot <- list()
    for (river in 1:length(HOBOS_dataset)) {
      mean_betw <- apply(ST_directed_betweennes_rivers[[river]],2,mean)
      sd_betw <- apply(ST_directed_betweennes_rivers[[river]],2,sd)
      
      n<- network(Simple_river_network[[river]], directed=T, diag=T)
      n %v% "B_values" <- mean_betw
      n %v% "Sd_values" <- mean_betw+sd_betw
      n %v% "edges_BC_values" <- mean_betw[1:length(mean_betw)-1]
      
      ST_betw_plot[[river]] <- ggplot(n, layout=as.matrix(Sites_coordinates[[river]][,3:4]),
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
    
    
    if(print.plots==TRUE){
      png(filename = paste(print.directory,"STbetweenness.png"),
          width = 715*6,height = 448*6,units = "px",res = 400)
      grid.arrange(arrangeGrob(grobs = ST_betw_plot,ncol=ceiling(length(HOBOS_dataset)/2)),
                   top="ST Dir NonW Betweenness")
      dev.off()}
    
  }
  
  ####_______________________________________________________________________
  # OUTPUTS _______________________####
  # Global matrix
  NonW_ST_matrix_rivers <- ST_matrix_rivers
  
  # Spatiotemporal matrix 
  NonW_ST_matrix_out_out <- ST_matrix_out_out
  NonW_ST_matrix_plots <- ST_matrix_plots
  
  # Spatiotemporal connectivity 
  NonW_ST_connectivity_value <- ST_connectivity_value
  NonW_ST_connectivity_plot <- ST_connectivity_plot
  
  NonW_ST_directed_Ocloseness_rivers <- list()
  NonW_ST_Oclo_plot <- list()
  NonW_ST_directed_Allcloseness_rivers <- list()
  NonW_ST_Allclo_plot <- list()
  NonW_ST_directed_betweennes_rivers <- list()
  NonW_ST_betw_plot <- list()
  
  if(Network_variables==T){
    # Spatiotemporal Out.closenness 
    NonW_ST_directed_Ocloseness_rivers <- ST_directed_Ocloseness_rivers
    NonW_ST_Oclo_plot <- ST_Oclo_plot
    
    # Spatiotemporal All.closenness 
    NonW_ST_directed_Allcloseness_rivers <- ST_directed_Allcloseness_rivers
    NonW_ST_Allclo_plot <- ST_Allclo_plot
    
    # Spatiotemporal Betweenness 
    NonW_ST_directed_betweennes_rivers <- ST_directed_betweennes_rivers
    NonW_ST_betw_plot <- ST_betw_plot 
  }
  
  Main_output <- list(Main_matrix=  NonW_ST_matrix_rivers,
                      STconmat=    NonW_ST_matrix_out_out,
                      STconmat_plo=NonW_ST_matrix_plots,
                      STcon=       NonW_ST_connectivity_value,
                      STcon_plo=   NonW_ST_connectivity_plot,
                      ST_Outcloseness=       NonW_ST_directed_Ocloseness_rivers,
                      ST_Outcloseness_plo=   NonW_ST_Oclo_plot,
                      ST_closeness=       NonW_ST_directed_Allcloseness_rivers,
                      ST_closeness_plo=   NonW_ST_Allclo_plot,
                      ST_Betweenness=       NonW_ST_directed_betweennes_rivers,
                      ST_Betweeness_plo=   NonW_ST_betw_plot)
  return(Main_output)
  ####_______________________________________________________________________
  
}# Function

