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
### 2 Old_HOBOS_calculation.R
### 3 SpaTemp_comparison.R
### 4 SpaTemp_Biolcomparison.R

# Nice colors? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")
setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 1. BIOLOGICAL Data uploading ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# IMPORTANT MESSAGE: Depending on which season you want to analye, you run different lines

# Spring dataset - May 2019
BDD <- read.csv2("BiolData/Matriz_primavera.csv", sep=";")
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

# Summer dataset - July 2019
BDD <- read.csv2("BiolData/Matriz_verano.csv", sep=";")
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

shan <- diversity(table_shannon%>%ungroup()%>%dplyr::select(-Code,-Riera),index = "shannon")

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

shan <- diversity(table_shannon%>%ungroup()%>%dplyr::select(-Code,-Riera),index = "shannon")

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
HOB_riv_ID<- Stream_order$ï..Riera
ups_dos <- Stream_order$UtoD

# Calculation of OLD HOBOS values to include them in the comparisons. 
source("Old_HOBOS_calculation.R")
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
Matching_HOBOS <- nn2(data=Sites_list_comb[,3:4], query = SampSites[,2:3], k=1, radius = 1)[[1]] 

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

names_scenarios <- c("Scenario A1",
                     "Scenario A2",
                     "Scenario B1",
                     "Scenario B2", 
                     "TotDur", "TotNum", "TotLeng")


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
    geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
    scale_fill_manual(values =  CunilleraPAL_corrected)+
    scale_color_manual(values =  CunilleraPAL_corrected)+
    xlab(paste("log(STcon)"))+
    ylab(paste("log(Richness)"))+
    labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
    labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
    theme_classic()+
    theme(legend.position="none")
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
      geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste("log(STcon)"))+
      ylab(paste("log(Shannon)"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none")
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
      geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste("log(STcon)"))+
      ylab(paste("log(Trait abundance)"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none")
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

# 5.4.1. Plot Bray ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
for (col_var in 1:ncol(STmatrix_BiolDissim%>%dplyr::select(-ID,
                                                    -BrayCurtis_Act, -Jaccard_Act,
                                                    -BrayCurtis_Pas, -Jaccard_Pas))){
  
  
  variable_x_temp <- STmatrix_BiolDissim%>%dplyr::select(-ID,
                                             -BrayCurtis_Act, -Jaccard_Act,
                                             -BrayCurtis_Pas, -Jaccard_Pas)

  variable_x_name <- colnames(variable_x_temp)[col_var]
  variable_x <- variable_x_temp[,col_var]
  
  variable_y_temp <- STmatrix_BiolDissim%>%dplyr::select(BrayCurtis_Act,BrayCurtis_Pas)
  plot_counter <- c(col_var, col_var+ncol(variable_x_temp))
  
  for (variable in 1:ncol(variable_y_temp)) {
    
    variable_y <- variable_y_temp[,variable]
    variable_y_model <- unlist(variable_y)
    variable_y_name <- colnames(variable_y_temp)
    
    factors <- STmatrix_BiolDissim%>%dplyr::select(ID)
    Id_random <- unlist(factors$ID)
    
    model <- lme(variable_y_model~log(variable_x+1), random = ~1|Id_random)
    results <- round(as.numeric(c(summary(model)[[4]]$fixed, summary(model)[[20]][2,5])),2)
    #Results table
    model_table <- round(summary(model)[[20]],3)
    rownames(model_table) <- c("Intercept", "STconmat")
    model_table<- rownames_to_column(as.data.frame(model_table))
    model_table<- cbind(Scenario=names_scenarios[col_var],Dispersion=Dispersal_group[variable],model_table)
    model_HOB_BDD_results[[plot_counter[variable]]] <- model_table 
    dataset <- cbind(factors, "X_var"=variable_x, variable_y)
    dataset$pred_values <- predict(model)   
    
    plots_HOB_BDD[[plot_counter[variable]]] <- ggplot(dataset,aes(x=log(X_var+1), y=variable_y,fill=ID,colour=ID))+
      geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
      geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, linetype=2)+
      geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste("log(STconmat)"))+
      ylab(paste("Bray curtis"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none")
    #facet_grid(.~ID)  
  }
}

model_result <- rbind(model_HOB_BDD_results[[1]],model_HOB_BDD_results[[2]],
                      model_HOB_BDD_results[[3]],model_HOB_BDD_results[[4]],
                      model_HOB_BDD_results[[5]],model_HOB_BDD_results[[6]],
                      model_HOB_BDD_results[[7]],model_HOB_BDD_results[[8]])
colnames(model_result)[[3]] <- "Fixed effects"
write.table(model_result, "Table_Results/Biol_BrayC_STconmat.txt", sep = ",")

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))

png(filename ="Figure/Biol_treat/BrayCurtis.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  legend_plots,
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=5, top="Bray Crutis")
dev.off()

plots_HOB_BDD_total[[4]] <- plots_HOB_BDD

# 5.4.2. Plot Jaccard ####
plots_HOB_BDD <- list()
model_HOB_BDD_results <- list()
for (col_var in 1:ncol(STmatrix_BiolDissim%>%dplyr::select(-ID,
                                                    -BrayCurtis_Act, -Jaccard_Act,
                                                    -BrayCurtis_Pas, -Jaccard_Pas))){
  
  
  variable_x_temp <- STmatrix_BiolDissim%>%dplyr::select(-ID,
                                                  -BrayCurtis_Act, -Jaccard_Act,
                                                  -BrayCurtis_Pas, -Jaccard_Pas)
  
  variable_x_name <- colnames(variable_x_temp)[col_var]
  variable_x <- variable_x_temp[,col_var]
  
  variable_y_temp <- STmatrix_BiolDissim%>%dplyr::select(Jaccard_Act,Jaccard_Pas)
  plot_counter <- c(col_var, col_var+ncol(variable_x_temp))
  
  Dispersal_group <- c("Active", "Passive")
  for (variable in 1:ncol(variable_y_temp)) {
    
    variable_y <- variable_y_temp[,variable]
    variable_y_model <- unlist(variable_y)
    variable_y_name <- colnames(variable_y_temp)
    
    factors <- STmatrix_BiolDissim%>%dplyr::select(ID)
    Id_random <- unlist(factors$ID)
    
    model <- lme(variable_y_model~log(variable_x+1), random = ~1|Id_random)
    results <- round(as.numeric(c(summary(model)[[4]]$fixed, summary(model)[[20]][2,5])),2)
    #Results table
    model_table <- round(summary(model)[[20]],3)
    rownames(model_table) <- c("Intercept", "STconmat")
    model_table<- rownames_to_column(as.data.frame(model_table))
    model_table<- cbind(Scenario=names_scenarios[col_var],Dispersion=Dispersal_group[variable],model_table)
    model_HOB_BDD_results[[plot_counter[variable]]] <- model_table 
    
    dataset <- cbind(factors, "X_var"=variable_x, variable_y)
    dataset$pred_values <- predict(model)   
    
    plots_HOB_BDD[[plot_counter[variable]]] <- ggplot(dataset,aes(x=log(X_var+1), y=variable_y,fill=ID,colour=ID))+
      geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
      geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, linetype=2)+
      geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
      scale_fill_manual(values =  CunilleraPAL_corrected)+
      scale_color_manual(values =  CunilleraPAL_corrected)+
      xlab(paste("log(STconmat)"))+
      ylab(paste("Jaccard"))+
      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
      labs(subtitle=paste(names_scenarios[col_var],"vs",variable_y_name[variable]))+
      theme_classic()+
      theme(legend.position="none")
    #facet_grid(.~ID)  
  }
}

model_result <- rbind(model_HOB_BDD_results[[1]],model_HOB_BDD_results[[2]],
                      model_HOB_BDD_results[[3]],model_HOB_BDD_results[[4]],
                      model_HOB_BDD_results[[5]],model_HOB_BDD_results[[6]],
                      model_HOB_BDD_results[[7]],model_HOB_BDD_results[[8]])
colnames(model_result)[[3]] <- "Fixed effects"
write.table(model_result, "Table_Results/Biol_Jacc_STconmat.txt", sep = ",")


legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=variable_y,fill=ID),shape=21, size=6)+
                            scale_fill_manual(values = CunilleraPAL_corrected, name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))

png(filename ="Figure/Biol_treat/Jaccard.png", 
    width = 1250*4, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  legend_plots,
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]], plots_HOB_BDD[[7]],plots_HOB_BDD[[8]],
  nrow=2,ncol=5, top="Jaccard")
dev.off()
plots_HOB_BDD_total[[5]] <- plots_HOB_BDD


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 6. FINAL RESULTS TABLES ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#### test chunk ####
#
#plots_HOB_BDD <- list()
#plots_HOB_BDD_ID_plot <- list()
#for (col_var in 1:ncol(HOB_BDD_match%>%select(-Riera, -Codi_HOBO,
#                                              -Latitud,-Longitud,
#                                              -ID, -DtoU,-ID_UpDo,-TotDur,-TotNum,-TotLeng,
#                                              -Samp_ID, -rich, -shan))){
#  ID_names <- unique(HOB_BDD_match$ID)
#  for (name_ID in 1:length(ID_names)) {
#    variable_x <- HOB_BDD_match%>%filter(ID==ID_names[name_ID])%>%select(-Riera, -Codi_HOBO,
#                                                                         -Latitud,-Longitud,
#                                                                         -ID, -DtoU,-ID_UpDo,-TotDur,-TotNum,-TotLeng,
#                                                                         -Samp_ID, -rich, -shan)
#    variable_x_name <- colnames(variable_x)[col_var]
#    variable_x <- variable_x[,col_var]
#    
#    variable_y <- HOB_BDD_match%>%filter(ID==ID_names[name_ID])%>%select(rich)
#    variable_y_model <- unlist(variable_y)
#    
#    factors <- HOB_BDD_match%>%filter(ID==ID_names[name_ID])%>%select(ID, DtoU,ID_UpDo)
#    Id_random <- unlist(factors$ID)
#    
#    model <- lm(log(variable_y_model+1)~log(variable_x+1))
#    if(is.na(model$coefficients[2])==TRUE){
#      results <- paste("-","-","-")
#    }else{
#      results <- round(as.numeric(c(summary(model)[[4]][,1], summary(model)[[4]][2,4])),2)
#    }
#    
#    col_fill <- CUNILLERA_pal("LGTBI")(length(ID_names))
#    
#    dataset <- cbind(factors, "X_var"=variable_x, variable_y)
#    #dataset$pred_values <- predict(model)
#    plots_HOB_BDD_ID_plot[[name_ID]]<- ggplot(dataset,aes(x=log(X_var+1), y=log(rich+1)))+
#      geom_point(color="grey30", alpha=0.5, shape=21, size=2,fill=col_fill[name_ID])+
#      geom_smooth(method = "lm", colour="grey60",alpha=0.1, se = TRUE, 
#                  linetype=2, fill=col_fill[name_ID])+
#      geom_smooth(method = "lm",alpha=0.2, se = T, fill="grey50", colour="black", size=2)+
#      xlab(paste("log(STcon)"))+
#      ylab(paste("log(richness)"))+
#      labs(caption = paste("Intercept=",results[1],"Slope=",results[2],"p-value=",results[3]))+
#      labs(title=paste(names_scenarios[col_var],"vs","richness"))+
#      theme_classic()+
#      theme(legend.position="none")
#  }
#  plots_HOB_BDD[[col_var]] <- plots_HOB_BDD_ID_plot
#}
#
#legend_plots<- get_legend(ggplot(dataset)+
#                            geom_point(aes(x=X_var,y=rich,fill=ID),shape=21, size=6)+
#                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
#                            theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical"))
#
#
#png(filename ="Figure/Biol_treat/Richness.png", 
#    width = 1000*4, height = 1000*4, 
#    units = "px",res = 300) 
#grid.arrange(
#  arrangeGrob(
#    plots_HOB_BDD[[1]][[1]],plots_HOB_BDD[[1]][[2]],plots_HOB_BDD[[1]][[3]],
#    plots_HOB_BDD[[1]][[4]],plots_HOB_BDD[[1]][[5]],plots_HOB_BDD[[1]][[6]], 
#    nrow=1,ncol=6),
#  
#  arrangeGrob(
#    plots_HOB_BDD[[2]][[1]],plots_HOB_BDD[[2]][[2]],plots_HOB_BDD[[2]][[3]],
#    plots_HOB_BDD[[2]][[4]],plots_HOB_BDD[[2]][[5]],plots_HOB_BDD[[2]][[6]],
#    nrow=1,ncol=6),
#  
#  arrangeGrob(
#    plots_HOB_BDD[[3]][[1]],plots_HOB_BDD[[3]][[2]],plots_HOB_BDD[[3]][[3]],
#    plots_HOB_BDD[[3]][[4]],plots_HOB_BDD[[3]][[5]],plots_HOB_BDD[[3]][[6]],
#    nrow=1,ncol=6),
#  
#  arrangeGrob(
#    plots_HOB_BDD[[4]][[1]],plots_HOB_BDD[[4]][[2]],plots_HOB_BDD[[4]][[3]],
#    plots_HOB_BDD[[4]][[4]],plots_HOB_BDD[[4]][[5]],plots_HOB_BDD[[4]][[6]],    
#    nrow=1,ncol=6),
#  arrangeGrob(legend_plots),
#  nrow=5,ncol=1)
#dev.off()