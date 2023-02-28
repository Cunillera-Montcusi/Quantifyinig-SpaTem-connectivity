# Load data management packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(nlme)
library(lsmeans)
library(multcomp)
library(car)
library(multcompView)
library(ade4)
library(vegan)
library(tidyverse)
detach("package:plyr", unload = TRUE)

### IMPORTANT: TO PROPERLY RUN THIS SCRIPT YOU NEED TO RUN the following scripts: 

## SpaTemp_HOBOS_treatment.R - The most important one as is the one calculating the ST indices and 
# selecting and creating the HOBOS dataset of interest and then preparing all the data in order to be analysed. 

## Old_HOBOS_calculation.R - Not that important but interesting. Here the Old HOBOS indices are calculated. This 
# indices are calcualted based on the amounts of "1" or "0" of the HOBOS table. They are interesting as they offer
# an extra perspective on the whole ST dynamics. 
## IMPORTANT: If you do not want to use them you will need to remove the parts where they are "called". Generally,
# they are named "Old_HOBOS" or similar names (TotDur, TotNum and TotLeng).

## SpaTemp_comparison.R - It is not essential to run this script. However, it is interesting to check how the 
# HOBOS obtained indices behave among them. Specially to understand what they mean and what they do not mean. So, 
# my recomendation is to run it. 

# SUGGESTED ORDER TO RUN the scripts: 
### 1 SpaTemp_HOBOS_treatment.R
### 2 FlowIntermittence_indices.R
### 3 SpaTemp_comparison.R
### 4 SpaTemp_Biolcomparison.R

# Nice colors? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")
setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1. BIOLOGICAL Data uploading ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Autumn dataset - November 2019
BDD <- read.csv2("BiolData/Matriz_otoño.csv", sep=";")
# We correct an error as there are two SA4 for October sampling, we sum them. 
BDD[which(BDD$Code=="SA4")[1],3:ncol(BDD)] <- apply(BDD[which(BDD$Code=="SA4"),3:ncol(BDD)],2,sum)
BDD <- BDD[-which(BDD$Code=="SA4")[2],]
BDD <- mutate(BDD, 
              Riera = case_when(
                str_detect(Code,"CA") ~ "C", 
                str_detect(Code,"H") ~ "VH",
                str_detect(Code,"MU") ~ "M",
                str_detect(Code,"R") ~ "R",
                str_detect(Code,"SA") ~ "SA",
                str_detect(Code,"SC") ~ "SC",
                str_detect(Code,"T") ~ "T",
                TRUE ~ "ERROR" )) 


Traits_val <- read.csv2("BiolData/traits.csv", row.names = 1)
bloc<-read.table("BiolData/trait.blocks.txt",h=T)

# 1.1 TRAITS implementation based on rel. abundances of each trait ####
# Filter species that are present or not present existing 
out <- c()
for (i in 1:length(rownames(Traits_val))) {
  temp_out <- which(rownames(Traits_val)==colnames(BDD%>%dplyr::select(-Date,-Riera,-Code,))[i])
  out[i] <- ifelse(length(temp_out)>0,temp_out,0)
}
Trait_matr_prep<-prep.fuzzy.var(Traits_val[out,],bloc$bloc)
apply(Trait_matr_prep,1,sum)

BDD_croped <- BDD%>%dplyr::select(matches(rownames(Trait_matr_prep)))
rownames(BDD_croped) <-BDD$Code

for(i in 1:nrow(BDD_croped)){BDD_croped[i,]<- BDD_croped[i,]/sum(BDD_croped[i,])}
apply(BDD_croped,1,sum)

Trait_abundances<-as.matrix(BDD_croped)%*%as.matrix(Trait_matr_prep)
  
Trait_abund_f4 <- data.frame("Code"=BDD$Code,"Date"=BDD$Date,"Riera"=BDD$Riera, "f4"=as.data.frame(Trait_abundances)$f4)
Trait_abund_f1 <- data.frame("Code"=BDD$Code,"Date"=BDD$Date,"Riera"=BDD$Riera, "f1"=as.data.frame(Trait_abundances)$f1)

# 1.2 TRAITS by filtering ####
# Traits addition by  filtering only species corresponding to the traits
Traits_val_edited <- data.frame("taxon"=row.names(Traits_val),Traits_val)

# Active dispersers -- f4==3 and 2
Trait_disp <- Traits_val_edited%>%dplyr::select(taxon,f4)%>%filter(f4==3)%>%
  rows_insert(Traits_val_edited%>%dplyr::select(taxon,f4)%>%filter(f4==2))
BDD_activ <- BDD%>%dplyr::select(Code, Date, Riera,matches(Trait_disp$taxon))

# Passive dispersers -- f1==3 and 2
Trait_disp <- Traits_val_edited%>%dplyr::select(taxon,f1)%>%filter(f1==3)%>%
  rows_insert(Traits_val_edited%>%dplyr::select(taxon,f1)%>%filter(f1==2))
BDD_pass <- BDD%>%dplyr::select(Code, Date, Riera,matches(Trait_disp$taxon))


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 2. DATA treatment and DIVERSITY values ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# 2.1 Active dispersers ####
Data_season <- BDD_activ %>%
                dplyr::select(-Date,-Riera)%>%
                pivot_longer(-Code)

colnames(Data_season) <- c("Code","Species","Abund")

#INDIVIDUAL diversity values active ________________________________________________________________________####
by_group_data_season <- group_by(Data_season, 
                                  Code, Species)

data_G_gamma <- by_group_data_season%>%
                group_by(Species)%>%
                summarise(n=sum(Abund))%>%
                mutate(s=1)%>%
                summarise(rich=sum(s))

data_A_alpha <- by_group_data_season%>%
                group_by(Code)%>%
                mutate(s=ifelse(Abund>0,1,0))%>%
                summarise(rich=sum(s))

library(vegan)
table_shannon <- by_group_data_season%>%
                 spread(key =Species,value = Abund, fill=0)%>%
                 bind_cols("Riera"=BDD_activ$Riera)

shan <- vegan::diversity(table_shannon%>%ungroup()%>%dplyr::select(-Code,-Riera),index = "shannon")

#SIMMILARITIES active ________________________________________________________________________####

table_shannon_PA <- by_group_data_season%>%
  mutate(PA=ifelse(Abund>0,1,0))%>%
  dplyr::select(-Abund)%>%
  spread(key =Species,value = PA, fill=0)%>%
  bind_cols("Riera"=BDD_activ$Riera)

Bray.distance <- list()
Jac.distance <- list()
for (riera in 1:length(unique(table_shannon$Riera))) {
  names_riera <- unique(table_shannon$Riera[order(table_shannon$Riera)])
  
  # Bray Curtis active
  Bray.distance[[riera]] <- vegdist(
                          decostand(table_shannon%>%
                                      filter(Riera==names_riera[riera])%>%ungroup()%>%dplyr::select(-Code,-Riera),
                          method = "hellinger"),
                          method = "bray")
  
  # Jaccard active
  Jac.distance[[riera]]<-vegdist(
                        decostand(table_shannon_PA%>%
                                    filter(Riera==names_riera[riera])%>%ungroup()%>%dplyr::select(-Code,-Riera),
                        method = "pa"),
                        method = "jaccard")
  }

names(Bray.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])
names(Jac.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])

BID_output_IND_act <- data.frame(Code=unique(by_group_data_season$Code),rich_f4=data_A_alpha$rich,shan_f4=shan,f4_abun=Trait_abund_f4$f4)
BID_output_DIST_act <- list(Bray.distance,Jac.distance)


# 2.2 Passive dispersers ####
Data_season <- BDD_pass %>%
  dplyr::select(-Date,-Riera)%>%
  pivot_longer(-Code)

colnames(Data_season) <- c("Code","Species","Abund")

#INDIVIDUAL diversity values passive ________________________________________________________________________####
by_group_data_season <- group_by(Data_season, 
                                 Code, Species)

data_G_gamma <- by_group_data_season%>%
  group_by(Species)%>%
  summarise(n=sum(Abund))%>%
  mutate(s=1)%>%
  summarise(rich=sum(s))

data_A_alpha <- by_group_data_season%>%
  group_by(Code)%>%
  mutate(s=ifelse(Abund>0,1,0))%>%
  summarise(rich=sum(s))

library(vegan)
table_shannon <- by_group_data_season%>%
  spread(key =Species,value = Abund, fill=0)%>%
  bind_cols("Riera"=BDD_pass$Riera)

shan <- vegan::diversity(table_shannon%>%ungroup()%>%dplyr::select(-Code,-Riera),index = "shannon")

#SIMMILARITIES passive ________________________________________________________________________####

table_shannon_PA <- by_group_data_season%>%
  mutate(PA=ifelse(Abund>0,1,0))%>%
  dplyr::select(-Abund)%>%
  spread(key =Species,value = PA, fill=0)%>%
  bind_cols("Riera"=BDD_pass$Riera)

Bray.distance <- list()
Jac.distance <- list()
for (riera in 1:length(unique(table_shannon$Riera))) {
  names_riera <- unique(table_shannon$Riera[order(table_shannon$Riera)])
  
  # Bray Curtis
  Bray.distance[[riera]] <- vegdist(
    decostand(table_shannon%>%
                filter(Riera==names_riera[riera])%>%ungroup()%>%dplyr::select(-Code,-Riera),
              method = "hellinger"),
    method = "bray")
  
  # Jaccard
  Jac.distance[[riera]]<-vegdist(
    decostand(table_shannon_PA%>%
                filter(Riera==names_riera[riera])%>%ungroup()%>%dplyr::select(-Code,-Riera),
              method = "pa"),
    method = "jaccard")
}

names(Bray.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])
names(Jac.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])

BID_output_IND_pas <- data.frame(Code=unique(by_group_data_season$Code),rich_f1=data_A_alpha$rich,shan_f1=shan,f1_abun=Trait_abund_f1$f1)
BID_output_DIST_pas <- list(Bray.distance,Jac.distance)


# BID values total
BID_output_IND <- left_join(BID_output_IND_act, BID_output_IND_pas, by="Code")

BID_output_DIST <- list("Active"=BID_output_DIST_act,"Passive"=BID_output_DIST_pas)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 3. CHARGING GEOGRAPHIC COORDINATES & MATCHING SAMPLES AND HOBOS ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# 3.1 We charge the field Samples and filter them according to the ones actually taken (from BDD)  ####
SampSites <- read.table("BiolData/Sites_1.txt",sep = ",",header = T)
out <- c()
for (match in 1:length(BDD$Code)) {
 a <- which(SampSites$Site==BDD$Code[match])
 out[match] <- ifelse(length(a)==0,0,a)
}
SampSites <- SampSites[out,]

# 3.2 We call the ST indices to MATCH them  ####
## You need to run the previous script in order to built the HOBOS dataset (a list with the HOBOS info for each river)
Sites_list # They are charged from SpaTemp_HOBOS_treatment.R
HOBOS_sites# They are charged from SpaTemp_HOBOS_treatment.R

# This two other objects are also built in the "SpaTemp_HOBOS_treatment.R" script. They contain the names 
## of each river and the order upstream to downstream within it. 
# They are used to Identify and locate river position between HOBOS and real Samples
HOB_riv_ID<- Stream_order$Riera
ups_dos <- Stream_order$UtoD

# Calculation of OLD HOBOS values to include them in the comparisons. 
source("FlowIntermittence_indices.R")
Old_HOBOS_comb

# 3.3 Geographical matching ####
# Converting Sampling SITES to apropriate geog. data 
str(SampSites)
a <- c()
b <- c()
for (posit in 1:nrow(SampSites)) {
a[posit] <-as.numeric(sub("(.{2})(.*)", "\\1.\\2", SampSites$lat[posit]))
b[posit] <-as.numeric(sub("(.{1})(.*)", "\\1.\\2", SampSites$long[posit]))
}
SampSites$lat<- a
SampSites$long<- b

# Converting HOBOS SITES to appropriate geog. data 
### *IMPORTANT STEP ####
#Here is where we include the calculated individual ST indices by "left_join"
Sites_list_comb <- rbind(Sites_list[[1]],Sites_list[[2]],Sites_list[[3]],Sites_list[[4]],Sites_list[[5]],
                         Sites_list[[6]],Sites_list[[7]] )
Sites_list_comb <- cbind(Sites_list_comb, "ID"=HOB_riv_ID,"DtoU"=ups_dos)%>%
                   mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))
Sites_list_comb <- Sites_list_comb%>%
                   left_join(NonW_ST_directed_out)%>%
                   left_join(WEIG_ST_directed_out)%>%
                   left_join(Un_NonW_ST_out)%>%
                   left_join(Un_WEIG_ST_out)%>%
                   left_join(Old_HOBOS_comb)
a <- c()
b <- c()
for (posit in 1:nrow(Sites_list_comb)) {
  a[posit] <-as.numeric(sub("(.{2})(.*)", "\\1.\\2", as.character(round(Sites_list_comb$Latitud[posit]))))
  b[posit] <-as.numeric(sub("(.{1})(.*)", "\\1.\\2", as.character(round(Sites_list_comb$Longitud[posit]))))
}
Sites_list_comb$Latitud<- a
Sites_list_comb$Longitud<- b

# Plot to check that the values match. 
par(mfrow=c(1,2))
plot(Sites_list_comb$Longitud, Sites_list_comb$Latitud)
plot(SampSites$long,SampSites$lat)
par(mfrow=c(1,1))

# Match the values
library(RANN) 
Matching_HOBOS <- nn2(data=Sites_list_comb[,4:3], query = SampSites[,3:2], k=1, radius = 1)[[1]] 

### *MANUAL CHECKING ####
# Here you must make sure that the assignation is correct between Samp_ID, ID, Riera, and DtoU
HOB_BDD_match <- as.data.frame(cbind(Sites_list_comb[Matching_HOBOS,], "Samp_ID"=SampSites$Site)) # Check assignation
data.frame(
HOB_BDD_match$Riera,
HOB_BDD_match$Codi_HOBO,
HOB_BDD_match$DtoU,
HOB_BDD_match$Samp_ID)

# Order the sites by upstream downstream directionality
HOB_BDD_match <- HOB_BDD_match%>%arrange(Riera,DtoU)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 4. MERGING ALL THE VALUES TOGETHER ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#Individual
HOB_BDD_match <- HOB_BDD_match%>%left_join(BID_output_IND, by=c("Samp_ID"="Code"))

#Matrices
### *IMPORTANT STEP ####
#Here is where we include the calculated matrix ST indices
NonW_Dir <- list()
WEIG_Dir<- list()
Un_NonW<- list()
Un_WEIG<- list()
for (riera in 1:length(unique(NonW_ST_matrix_out_out))) {
  names_riera <- unique(Sites_list_comb$ID)
  
  values_places <- unlist(HOB_BDD_match%>%filter(ID==names_riera[riera])%>%dplyr::select(DtoU))
  print(values_places)
  
  NonW_Dir[[riera]] <- as.dist(t(NonW_ST_matrix_out_out[[riera]][values_places,values_places]))
  WEIG_Dir[[riera]] <- as.dist(t(WEIG_ST_matrix_out_out[[riera]][values_places,values_places]))
  Un_NonW[[riera]] <- as.dist(t(Un_NonW_ST_matrix_out_out[[riera]][values_places,values_places]))
  Un_WEIG[[riera]] <- as.dist(t(Un_WEIG_ST_matrix_out_out[[riera]][values_places,values_places]))
}

### *MANUAL CHECKING ####
# We eliminate R and G matrices due to the lack of sites for this two sampling places
## For diversities (Bray and Jaccard) eliminate positions 2 and 4
## For ST matrices eliminate positions 3
HOB_BDD_match_matrix <- list("BrayCurtis_Act"=BID_output_DIST$Active[[1]][-c(2,4)],
                             "Jaccard_Act"=BID_output_DIST$Active[[2]][-c(2,4)],
                             "BrayCurtis_Pas"=BID_output_DIST$Passive[[1]][-c(2,4)],
                             "Jaccard_Pas"=BID_output_DIST$Passive[[2]][-c(2,4)],
                             "NonW_Dir"=NonW_Dir[-c(3)], 
                             "WEIG_Dir"=WEIG_Dir[-c(3)], 
                             "Un_NonW"=Un_NonW[-c(3)], 
                             "Un_WEIG"=Un_WEIG[-c(3)])

# We eliminate R and G in Individual due to the lack of samples from this period
HOB_BDD_match <- HOB_BDD_match%>%filter(Samp_ID !=c("R4","G1"))

# FINAL MERGED VALUES - INDIVIDUAL AND DISTANCE (MATRIX) BASED 
HOB_BDD_match
HOB_BDD_match_matrix

### *SMALL EDITIONS TO ELIMINATE Ocl/Acl/Bet ####
HOB_BDD_match <- HOB_BDD_match%>%dplyr::select(-c(NonW_Dir_Ocl,NonW_Dir_Acl,NonW_Dir_Bet,
                                           WEIG_Dir_Ocl,WEIG_Dir_Acl,WEIG_Dir_Bet,
                                           Un_NonW_Ocl,Un_NonW_Acl,Un_NonW_Bet,
                                           Un_WEIG_Ocl,Un_WEIG_Acl,Un_WEIG_Bet))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 5. Relationships with BIOL DATA ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

CunilleraPAL_corrected <- CUNILLERA_pal("LGTBI")(7)[-3]

names_scenarios <- c("DirBin",
                     "DirWei",
                     "UndBin",
                     "UndWei", 
                     "TotDur", "TotNum", "TotLeng")
axistextcolours <- c(viridis::viridis(5)[3],viridis::inferno(5)[3],
                     viridis::viridis(5)[3],viridis::inferno(5)[3])

axisXtitles <- c(c("Isolated   <-- log(STcon) -->  Connected"), 
                 c("Low disp. resist <-- log(STcon) --> High disp. resist"),
                 c("Isolated   <-- log(STcon) -->  Connected"), 
                 c("Low disp. resist <-- log(STcon) --> High disp. resist"))

# 5.1. Richness ####
plots_HOB_BDD_total <- list()
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
for (col_var in 1:ncol(HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                              -Latitud,-Longitud,
                                              -ID, -DtoU,-ID_UpDo,
                                              -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                              -rich_f4, -shan_f4,-f4_abun,
                                              -rich_f1,-shan_f1,-f1_abun))){

  
  variable_x_temp <- HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                       -Latitud,-Longitud,
                                       -ID, -DtoU,-ID_UpDo,
                                       -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                       -rich_f4, -shan_f4,-f4_abun,
                                       -rich_f1,-shan_f1,-f1_abun)

  variable_x_name <- colnames(variable_x_temp)[col_var]
  variable_x <- variable_x_temp[,col_var]
  
  variable_y_temp <- HOB_BDD_match%>%dplyr::select(rich_f4,rich_f1)
  plot_counter <- c(col_var, col_var+ncol(variable_x_temp))
  
for (variable in 1:ncol(variable_y_temp)) {
    
  variable_y <- variable_y_temp[,variable]
  variable_y_model <- unlist(variable_y)
  variable_y_name <- colnames(variable_y_temp)
  
  factors <- HOB_BDD_match%>%dplyr::select(ID, DtoU,ID_UpDo)
  Id_random <- unlist(factors$ID)
  
  model <- lme(log(variable_y_model+1)~log(variable_x+1), random = ~1|Id_random)
  results <- round(as.numeric(c(summary(model)[[4]]$fixed, summary(model)[[20]][2,5])),2)
  #Results table
  model_table <- round(summary(model)[[20]],3)
  rownames(model_table) <- c("Intercept", "STconmat")
  model_table<- rownames_to_column(as.data.frame(model_table))
  model_table<- cbind(Scenario=names_scenarios[col_var],Dispersion=Dispersal_group[variable],model_table)
  model_HOB_BDD_results[[plot_counter[variable]]] <- model_table 
  
  dataset <- cbind(factors, "X_var"=variable_x, variable_y)
  dataset$pred_values <- predict(model)   
  
plots_HOB_BDD[[plot_counter[variable]]] <- ggplot(dataset,aes(x=log(X_var+1), y=log(variable_y+1),fill=ID,colour=ID))+
    geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
    geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, linetype=2)+
    geom_abline(slope =results[2],intercept = results[1], colour="black", size=2)+
    #geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
    scale_fill_manual(values =  CunilleraPAL_corrected)+
    scale_color_manual(values =  CunilleraPAL_corrected)+
    xlab(paste(axisXtitles[col_var]))+
    ylab(paste("log(Richness)"))+
    labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
    labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
    theme_classic()+
    theme(legend.position="none",
          axis.line.x.bottom=element_line(color=axistextcolours[col_var], size=3),
          axis.text.x = element_text(color =axistextcolours[col_var] ),
          axis.title.x = element_text(color =axistextcolours[col_var] , size=8 ))
    #facet_grid(.~ID)  
  }
}

model_result <- rbind(model_HOB_BDD_results[[1]],model_HOB_BDD_results[[2]],
                      model_HOB_BDD_results[[3]],model_HOB_BDD_results[[4]],
                      model_HOB_BDD_results[[5]],model_HOB_BDD_results[[6]],
                      model_HOB_BDD_results[[7]],model_HOB_BDD_results[[8]])
colnames(model_result)[[3]] <- "Fixed effects"
write.table(model_result, "Table_Results/Biol_Rich_STcon.txt", sep = ",")

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y ,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))


png(filename ="Figure/Biol_treat/Richness.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  legend_plots,
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=5, top="Richness")
dev.off()
plots_HOB_BDD_total[[1]] <- plots_HOB_BDD


# 5.2. Shannon ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
for (col_var in 1:ncol(HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                              -Latitud,-Longitud,
                                              -ID, -DtoU,-ID_UpDo,
                                              -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                              -rich_f4, -shan_f4,-f4_abun,
                                              -rich_f1,-shan_f1,-f1_abun))){
  
  
  variable_x_temp <- HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                            -Latitud,-Longitud,
                                            -ID, -DtoU,-ID_UpDo,
                                            -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                            -rich_f4, -shan_f4,-f4_abun,
                                            -rich_f1,-shan_f1,-f1_abun)
  
  variable_x_name <- colnames(variable_x_temp)[col_var]
  variable_x <- variable_x_temp[,col_var]
  
  variable_y_temp <- HOB_BDD_match%>%dplyr::select(shan_f4,shan_f1)
  plot_counter <- c(col_var, col_var+ncol(variable_x_temp))
  
  for (variable in 1:ncol(variable_y_temp)) {
    
    variable_y <- variable_y_temp[,variable]
    variable_y_model <- unlist(variable_y)
    variable_y_name <- colnames(variable_y_temp)
    
    factors <- HOB_BDD_match%>%dplyr::select(ID, DtoU,ID_UpDo)
    Id_random <- unlist(factors$ID)
    
    model <- lme(log(variable_y_model+1)~log(variable_x+1), random = ~1|Id_random)
    results <- round(as.numeric(c(summary(model)[[4]]$fixed, summary(model)[[20]][2,5])),2)
    #Results table
    model_table <- round(summary(model)[[20]],3)
    rownames(model_table) <- c("Intercept", "STconmat")
    model_table<- rownames_to_column(as.data.frame(model_table))
    model_table<- cbind(Scenario=names_scenarios[col_var],Dispersion=Dispersal_group[variable],model_table)
    model_HOB_BDD_results[[plot_counter[variable]]] <- model_table 
    
    dataset <- cbind(factors, "X_var"=variable_x, variable_y)
    dataset$pred_values <- predict(model)   
    
    plots_HOB_BDD[[plot_counter[variable]]] <- ggplot(dataset,aes(x=log(X_var+1), y=log(variable_y+1),fill=ID,colour=ID))+
      geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
      geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, linetype=2)+
      geom_abline(slope =results[2],intercept = results[1], colour="black", size=2)+
      #geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste(axisXtitles[col_var]))+
      ylab(paste("log(Shannon)"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none",
            axis.line.x.bottom=element_line(color=axistextcolours[col_var], size=3),
            axis.text.x = element_text(color =axistextcolours[col_var] ),
            axis.title.x = element_text(color =axistextcolours[col_var], size=8  ))
    #facet_grid(.~ID)  
  }
}

model_result <- rbind(model_HOB_BDD_results[[1]],model_HOB_BDD_results[[2]],
                      model_HOB_BDD_results[[3]],model_HOB_BDD_results[[4]],
                      model_HOB_BDD_results[[5]],model_HOB_BDD_results[[6]],
                      model_HOB_BDD_results[[7]],model_HOB_BDD_results[[8]])
colnames(model_result)[[3]] <- "Fixed effects"
write.table(model_result, "Table_Results/Biol_ShaWei_STcon.txt", sep = ",")

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y ,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))


png(filename ="Figure/Biol_treat/Shannon.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  legend_plots,
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=5, top="Shannon")
dev.off()

plots_HOB_BDD_total[[2]] <- plots_HOB_BDD


# 5.3. Trait abund. ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
for (col_var in 1:ncol(HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                              -Latitud,-Longitud,
                                              -ID, -DtoU,-ID_UpDo,
                                              -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                              -rich_f4, -shan_f4,-f4_abun,
                                              -rich_f1,-shan_f1,-f1_abun))){
  
  
  variable_x_temp <- HOB_BDD_match%>%dplyr::select(-Riera, -Codi_HOBO,
                                            -Latitud,-Longitud,
                                            -ID, -DtoU,-ID_UpDo,
                                            -Samp_ID, -TotDur, -TotNum, -TotLeng, 
                                            -rich_f4, -shan_f4,-f4_abun,
                                            -rich_f1,-shan_f1,-f1_abun)
  
  variable_x_name <- colnames(variable_x_temp)[col_var]
  variable_x <- variable_x_temp[,col_var]
  
  variable_y_temp <- HOB_BDD_match%>%dplyr::select(f4_abun,f1_abun)
  plot_counter <- c(col_var, col_var+ncol(variable_x_temp))
  
  for (variable in 1:ncol(variable_y_temp)) {
    
    variable_y <- variable_y_temp[,variable]
    variable_y_model <- unlist(variable_y)
    variable_y_name <- colnames(variable_y_temp)
    
    factors <- HOB_BDD_match%>%dplyr::select(ID, DtoU,ID_UpDo)
    Id_random <- unlist(factors$ID)
    
    model <- lme(log(variable_y_model+1)~log(variable_x+1), random = ~1|Id_random)
    results <- round(as.numeric(c(summary(model)[[4]]$fixed, summary(model)[[20]][2,5])),2)
    #Results table
    model_table <- round(summary(model)[[20]],3)
    rownames(model_table) <- c("Intercept", "STconmat")
    model_table<- rownames_to_column(as.data.frame(model_table))
    model_table<- cbind(Scenario=names_scenarios[col_var],Dispersion=Dispersal_group[variable],model_table)
    model_HOB_BDD_results[[plot_counter[variable]]] <- model_table 
    
    dataset <- cbind(factors, "X_var"=variable_x, variable_y)
    dataset$pred_values <- predict(model)   
    
    plots_HOB_BDD[[plot_counter[variable]]] <- ggplot(dataset,aes(x=log(X_var+1), y=log(variable_y+1),fill=ID,colour=ID))+
      geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
      geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, linetype=2)+
      geom_abline(slope =results[2],intercept = results[1], colour="black", size=2)+
      #geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste(axisXtitles[col_var]))+
      ylab(paste("log(Trait abundance)"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none",
            axis.line.x.bottom=element_line(color=axistextcolours[col_var], size=3),
            axis.text.x = element_text(color =axistextcolours[col_var] ),
            axis.title.x = element_text(color =axistextcolours[col_var], size=8  ))
    #facet_grid(.~ID)  
  }
}

model_result <- rbind(model_HOB_BDD_results[[1]],model_HOB_BDD_results[[2]],
                      model_HOB_BDD_results[[3]],model_HOB_BDD_results[[4]],
                      model_HOB_BDD_results[[5]],model_HOB_BDD_results[[6]],
                      model_HOB_BDD_results[[7]],model_HOB_BDD_results[[8]])
colnames(model_result)[[3]] <- "Fixed effects"
write.table(model_result, "Table_Results/Biol_TraitAb_STcon.txt", sep = ",")

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y ,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))


png(filename ="Figure/Biol_treat/TraitAbun.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  legend_plots,
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=5, top="Trait abundance")
dev.off()

plots_HOB_BDD_total[[3]] <- plots_HOB_BDD


# 5.4. Matching with Bray & Jaccard ####
output <- matrix(nrow=length(unlist(HOB_BDD_match_matrix$BrayCurtis_Act)))
for (lis_elem in 1:length(HOB_BDD_match_matrix)) {
  out <- c()
  for (river in 1:length(HOB_BDD_match_matrix$BrayCurtis_Act)) {
    col_names_river <- names(HOB_BDD_match_matrix$BrayCurtis_Act)
    
    out_temp <- c()
    out_temp <- cbind(as.matrix(HOB_BDD_match_matrix[[lis_elem]][[river]])[
      upper.tri(as.matrix(HOB_BDD_match_matrix[[lis_elem]][[river]]))])
    out_temp <- as.data.frame(cbind("ID"=rep(col_names_river[river],nrow(out_temp)),out_temp))
    out<- rbind(out, out_temp)
  }
  colnames(out)[2] <- names(HOB_BDD_match_matrix)[lis_elem]
  output[,1] <- out[,1]
  output <- as.data.frame(cbind(output, out))
  output <- output%>%dplyr::select(-ID)
  output[,lis_elem+1] <- as.numeric(output[,lis_elem+1])
}
colnames(output)[1] <- "ID"

STmatrix_BiolDissim <- output

axisXtitles <- c(c("Isolated   <-- STconmat -->  Connected"), 
                 c("Low disp. resist <-- STconmat --> High disp. resist"),
                 c("Isolated   <-- STconmat -->  Connected"), 
                 c("Low disp. resist <-- STconmat --> High disp. resist"))

leng_dist_matrix <- ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[1]]))+
ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[2]]))+
ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[3]]))+
ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[4]]))+
ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[5]]))+
ncol(as.matrix(HOB_BDD_match_matrix$BrayCurtis_Act[[6]]))

matrices_total_analyse <- list()
for (matrices_tot in 1:length(HOB_BDD_match_matrix)) {
Matrix_Beta <- HOB_BDD_match_matrix[[matrices_tot]]
full_matrix <- matrix(nrow = leng_dist_matrix,ncol = leng_dist_matrix,data = 0)
for (matrix_number in 1:length(Matrix_Beta)) {
ncol_leng <- ncol(as.matrix(Matrix_Beta[[matrix_number]]))
mac_ncol_leng <- which(apply(full_matrix,2,sum)==0)[ncol_leng]
min_ncol_leng <- which(apply(full_matrix,2,sum)==0)[1]
full_matrix[min_ncol_leng:mac_ncol_leng,
            min_ncol_leng:mac_ncol_leng] <- as.matrix(Matrix_Beta[[matrix_number]])
}
zeros <- which(full_matrix==0)
full_matrix[zeros] <- mean(full_matrix[-zeros])
matrices_total_analyse[[matrices_tot]] <- full_matrix
}
names(matrices_total_analyse) <- names(HOB_BDD_match_matrix)


# 5.4.1. Plot Bray ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
variable_y_name <- c("BrayCurtis_Act","BrayCurtis_Pas")
Dispersal_group <- c("Active", "Passive")
# BrayCurtis_Act
out_results <- data.frame()
for (stconmat_val in 1:4) {
a <- vegan::mantel(ydis = as.dist(matrices_total_analyse$BrayCurtis_Act),
                   xdis = as.dist(matrices_total_analyse[[4+stconmat_val]]),permutations = 1000,method = "spearman")
results <- data.frame("Group"=Dispersal_group[1],
                      "Scenario"=names_scenarios[stconmat_val],
                      "Mantel coefficient"=round(a$statistic,3), "p.value"=round(a$signif,3))

plots_HOB_BDD[[stconmat_val]] <- data.frame(
"X"=matrices_total_analyse[[4+stconmat_val]][lower.tri(matrices_total_analyse[[4+stconmat_val]])],
"Y"=matrices_total_analyse$BrayCurtis_Act[lower.tri(matrices_total_analyse$BrayCurtis_Act)]) %>% 
  ggplot(aes(x=X,y=Y))+
  geom_point(color="grey20",fill="grey40",alpha=0.8, shape=21, size=2)+
  geom_smooth(method = "lm", colour="black",alpha=0.1, se = TRUE, linetype=1, linewidth=2)+
  xlab(paste(axisXtitles[stconmat_val]))+
  ylab(paste("Bray curtis"))+
  labs(caption = paste("Mantel coefficient=",results[3],"p.value=",results[4]))+
  labs(subtitle=paste(names_scenarios[stconmat_val],"vs",variable_y_name[1]))+
  theme_classic()+
  theme(legend.position="none",
        axis.line.x.bottom=element_line(color=axistextcolours[stconmat_val], linewidth=3),
        axis.text.x = element_text(color =axistextcolours[stconmat_val] ),
        axis.title.x = element_text(color =axistextcolours[stconmat_val] , size=8 ))
out_results <- bind_rows(out_results, results)
}

# BrayCurtis_Pas
for (stconmat_val in 1:4) {
a <- vegan::mantel(ydis = as.dist(matrices_total_analyse$BrayCurtis_Pas),
                   xdis = as.dist((matrices_total_analyse[[4+stconmat_val]])),permutations = 1000,method = "spearman")
results <- data.frame("Group"=Dispersal_group[2],
                      "Scenario"=names_scenarios[stconmat_val],
                      "Mantel coefficient"=round(a$statistic,3), "p.value"=round(a$signif,3))

plots_HOB_BDD[[stconmat_val+4]] <- data.frame(
  "X"=(matrices_total_analyse[[4+stconmat_val]][lower.tri(matrices_total_analyse[[4+stconmat_val]])]),
  "Y"=matrices_total_analyse$BrayCurtis_Pas[lower.tri(matrices_total_analyse$BrayCurtis_Pas)]) %>% 
  ggplot(aes(x=X,y=Y))+
  geom_point(color="grey20",fill="grey40",alpha=0.8, shape=21, size=2)+
  geom_smooth(method = "lm", colour="black",alpha=0.1, se = TRUE, linetype=1, linewidth=2)+
  xlab(paste(axisXtitles[stconmat_val]))+
  ylab(paste("Bray curtis"))+
  labs(caption = paste("Mantel coefficient=",results[3],"p.value=",results[4]))+
  labs(subtitle=paste(names_scenarios[stconmat_val],"vs",variable_y_name[2]))+
  theme_classic()+
  theme(legend.position="none",
        axis.line.x.bottom=element_line(color=axistextcolours[stconmat_val], linewidth=3),
        axis.text.x = element_text(color =axistextcolours[stconmat_val] ),
        axis.title.x = element_text(color =axistextcolours[stconmat_val] , size=8 ))
out_results <- bind_rows(out_results, results)
}

model_result <- out_results
write.table(model_result, "Table_Results/Biol_BrayC_STconmat.txt", sep = ",")

png(filename ="Figure/Biol_treat/BrayCurtis.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=4, top="Bray Crutis")
dev.off()
plots_HOB_BDD_total[[4]] <- plots_HOB_BDD

# 5.4.2. Plot Jaccard ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
variable_y_name <- c("Jaccard_Act","Jaccard_Pas")
Dispersal_group <- c("Active", "Passive")
# Jaccard_Act
out_results <- data.frame()
for (stconmat_val in 1:4) {
  a <- vegan::mantel(ydis = as.dist(matrices_total_analyse$Jaccard_Act),
                     xdis = as.dist(matrices_total_analyse[[4+stconmat_val]]),permutations = 1000,method = "spearman")
  results <- data.frame("Group"=Dispersal_group[1],
                        "Scenario"=names_scenarios[stconmat_val],
                        "Mantel coefficient"=round(a$statistic,3), "p.value"=round(a$signif,3))
  
  plots_HOB_BDD[[stconmat_val]] <- data.frame(
    "X"=matrices_total_analyse[[4+stconmat_val]][lower.tri(matrices_total_analyse[[4+stconmat_val]])],
    "Y"=matrices_total_analyse$Jaccard_Act[lower.tri(matrices_total_analyse$Jaccard_Act)]) %>% 
    ggplot(aes(x=X,y=Y))+
    geom_point(color="grey20",fill="grey40",alpha=0.8, shape=21, size=2)+
    geom_smooth(method = "lm", colour="black",alpha=0.1, se = TRUE, linetype=1, linewidth=2)+
    xlab(paste(axisXtitles[stconmat_val]))+
    ylab(paste("Jaccard"))+
    labs(caption = paste("Mantel coefficient=",results[3],"p.value=",results[4]))+
    labs(subtitle=paste(names_scenarios[stconmat_val],"vs",variable_y_name[1]))+
    theme_classic()+
    theme(legend.position="none",
          axis.line.x.bottom=element_line(color=axistextcolours[stconmat_val], linewidth=3),
          axis.text.x = element_text(color =axistextcolours[stconmat_val] ),
          axis.title.x = element_text(color =axistextcolours[stconmat_val] , size=8 ))
  out_results <- bind_rows(out_results, results)
}

# Jaccard_Pas
for (stconmat_val in 1:4) {
  a <- vegan::mantel(ydis = as.dist(matrices_total_analyse$Jaccard_Pas),
                     xdis = as.dist((matrices_total_analyse[[4+stconmat_val]])),permutations = 1000,method = "spearman")
  results <- data.frame("Group"=Dispersal_group[2],
                        "Scenario"=names_scenarios[stconmat_val],
                        "Mantel coefficient"=round(a$statistic,3), "p.value"=round(a$signif,3))
  
  plots_HOB_BDD[[stconmat_val+4]] <- data.frame(
    "X"=(matrices_total_analyse[[4+stconmat_val]][lower.tri(matrices_total_analyse[[4+stconmat_val]])]),
    "Y"=matrices_total_analyse$Jaccard_Pas[lower.tri(matrices_total_analyse$Jaccard_Pas)]) %>% 
    ggplot(aes(x=X,y=Y))+
    geom_point(color="grey20",fill="grey40",alpha=0.8, shape=21, size=2)+
    geom_smooth(method = "lm", colour="black",alpha=0.1, se = TRUE, linetype=1, linewidth=2)+
    xlab(paste(axisXtitles[stconmat_val]))+
    ylab(paste("Jaccard"))+
    labs(caption = paste("Mantel coefficient=",results[3],"p.value=",results[4]))+
    labs(subtitle=paste(names_scenarios[stconmat_val],"vs",variable_y_name[2]))+
    theme_classic()+
    theme(legend.position="none",
          axis.line.x.bottom=element_line(color=axistextcolours[stconmat_val], linewidth=3),
          axis.text.x = element_text(color =axistextcolours[stconmat_val] ),
          axis.title.x = element_text(color =axistextcolours[stconmat_val] , size=8 ))
  out_results <- bind_rows(out_results, results)
}

model_result <- out_results
write.table(model_result, "Table_Results/Biol_Jacc_STconmat.txt", sep = ",")

png(filename ="Figure/Biol_treat/Jaccard.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=4, top="Jaccard")
dev.off()
plots_HOB_BDD_total[[5]] <- plots_HOB_BDD


