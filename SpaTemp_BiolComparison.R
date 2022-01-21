# Load data management packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(nlme)
library(lsmeans)
library(multcomp)
library(car)
library(multcompView)

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
### 4 SpaTemp_comparison.R

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


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 2. DATA treatment and DIVERSITY values ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
Data_season <- BDD %>%
                select(-Date,-Riera)%>%
                pivot_longer(-Code)

colnames(Data_season) <- c("Code","Species","Abund")

#RICHNESS ________________________________________________________________________####
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
                 bind_cols("Riera"=BDD$Riera)

shan <- diversity(table_shannon%>%ungroup()%>%select(-Code,-Riera),index = "shannon")

#SIMMILARITIES ________________________________________________________________________####

table_shannon_PA <- by_group_data_season%>%
  mutate(PA=ifelse(Abund>0,1,0))%>%
  select(-Abund)%>%
  spread(key =Species,value = PA, fill=0)%>%
  bind_cols("Riera"=BDD$Riera)

Bray.distance <- list()
Jac.distance <- list()
for (riera in 1:length(unique(table_shannon$Riera))) {
  names_riera <- unique(table_shannon$Riera[order(table_shannon$Riera)])
  
  # Bray Curtis
  Bray.distance[[riera]] <- vegdist(
                          decostand(table_shannon%>%
                                      filter(Riera==names_riera[riera])%>%ungroup()%>%select(-Code,-Riera),
                          method = "hellinger"),
                          method = "bray")
  
  # Jaccard
  Jac.distance[[riera]]<-vegdist(
                        decostand(table_shannon_PA%>%
                                    filter(Riera==names_riera[riera])%>%ungroup()%>%select(-Code,-Riera),
                        method = "pa"),
                        method = "jaccard")
  }

names(Bray.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])
names(Jac.distance) <- unique(table_shannon$Riera[order(table_shannon$Riera)])

BID_output_IND <- data.frame(Code=unique(by_group_data_season$Code),rich=data_A_alpha$rich,shan)

BID_output_DIST <- list(Bray.distance,Jac.distance)

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
HOB_riv_ID
ups_dos

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
str(Sites_list_comb)
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
for (riera in 1:length(unique(HOB_BDD_match$ID))) {
  names_riera <- unique(HOB_BDD_match$ID)
  
  values_places <- unlist(HOB_BDD_match%>%filter(ID==names_riera[riera])%>%select(DtoU))
  
  NonW_Dir[[riera]] <- as.dist(t(NonW_ST_matrix_out_out[[riera]][values_places,values_places]))
  WEIG_Dir[[riera]] <- as.dist(t(WEIG_ST_matrix_out_out[[riera]][values_places,values_places]))
  Un_NonW[[riera]] <- as.dist(t(Un_NonW_ST_matrix_out_out[[riera]][values_places,values_places]))
  Un_WEIG[[riera]] <- as.dist(t(Un_WEIG_ST_matrix_out_out[[riera]][values_places,values_places]))
}

### *MANUAL CHECKING ####
# We eliminate R and G matrices due to the lack of sites for this two sampling places
## For diversities (Bray and Jaccard) eliminate positions 2 and 4
## For ST matrices eliminate positions 3
HOB_BDD_match_matrix <- list("BrayCurtis"=BID_output_DIST[[1]][-c(2,4)],
                             "Jaccard"=BID_output_DIST[[2]][-c(2,4)],
                             "NonW_Dir"=NonW_Dir[-c(3)], 
                             "WEIG_Dir"=WEIG_Dir[-c(3)], 
                             "Un_NonW"=Un_NonW[-c(3)], 
                             "Un_WEIG"=Un_WEIG[-c(3)])

# We eliminate R and G in Individual due to the lack of samples from this period
HOB_BDD_match <- HOB_BDD_match%>%filter(Samp_ID !=c("R4","G1"))

# FINAL MERGED VALUES - INDIVIDUAL AND DISTANCE (MATRIX) BASED 
HOB_BDD_match
HOB_BDD_match_matrix

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 5. Relationships with BIOL DATA ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# 5.1. Richness ####
plots_HOB_BDD <- list()
for (col_var in 1:ncol(HOB_BDD_match%>%select(-Riera, -Codi_HOBO,
                                              -Latitud,-Longitud,
                                              -ID, -DtoU,-ID_UpDo,
                                              -Samp_ID, -rich, -shan))){
  
variable_x <- HOB_BDD_match%>%select(-Riera, -Codi_HOBO,
                                     -Latitud,-Longitud,
                                     -ID, -DtoU,-ID_UpDo,
                                     -Samp_ID, -rich, -shan)
variable_x_name <- colnames(variable_x)[col_var]
variable_x <- variable_x[,col_var]

variable_y <- HOB_BDD_match%>%select(rich)

factors <- HOB_BDD_match%>%select(ID, DtoU,ID_UpDo)
  
dataset <- cbind(factors, "X_var"=variable_x, variable_y)
plots_HOB_BDD[[col_var]] <- ggplot(dataset,aes(x=X_var, y=rich,fill=ID))+
                            geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
                            geom_smooth(aes(colour=ID),method = "lm",alpha=0.1, se = F, linetype=2)+
                            geom_smooth(method = "lm",alpha=0.1, se = T, fill="grey50", colour="black", size=2)+
                            scale_fill_CUNILLERA(palette = "LGTBI")+
                            scale_color_CUNILLERA(palette = "LGTBI")+
                            xlab(variable_x_name)+
                            labs(title=paste(variable_x_name,"vs","Richness"))+
                            theme_classic()+
                            theme(legend.position="none")
}

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=rich,fill=ID),shape=21, size=6)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical"))


png(filename ="Figure/Biol_treat/Richness.png", 
    width = 1000*4, height = 1000*4, 
    units = "px",res = 300) 
grid.arrange(
plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]], 
plots_HOB_BDD[[5]],plots_HOB_BDD[[6]],plots_HOB_BDD[[7]],plots_HOB_BDD[[8]], 
plots_HOB_BDD[[9]],plots_HOB_BDD[[10]],plots_HOB_BDD[[11]],plots_HOB_BDD[[12]], 
plots_HOB_BDD[[13]],plots_HOB_BDD[[14]],plots_HOB_BDD[[15]],plots_HOB_BDD[[16]], 
plots_HOB_BDD[[17]],plots_HOB_BDD[[18]],plots_HOB_BDD[[19]], legend_plots,
nrow=5,ncol=4)
dev.off()

# 5.2. Shannon ####
plots_HOB_BDD <- list()
for (col_var in 1:ncol(HOB_BDD_match%>%select(-Riera, -Codi_HOBO,
                                              -Latitud,-Longitud,
                                              -ID, -DtoU,-ID_UpDo,
                                              -Samp_ID, -rich, -shan))){
  
  variable_x <- HOB_BDD_match%>%select(-Riera, -Codi_HOBO,
                                       -Latitud,-Longitud,
                                       -ID, -DtoU,-ID_UpDo,
                                       -Samp_ID, -rich, -shan)
  variable_x_name <- colnames(variable_x)[col_var]
  variable_x <- variable_x[,col_var]
  
  variable_y <- HOB_BDD_match%>%select(shan)
  
  factors <- HOB_BDD_match%>%select(ID, DtoU,ID_UpDo)
  
  dataset <- cbind(factors, "X_var"=variable_x, variable_y)
  plots_HOB_BDD[[col_var]] <- ggplot(dataset,aes(x=X_var, y=shan,fill=ID))+
    geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
    geom_smooth(aes(colour=ID),method = "lm",alpha=0.1, se = F, linetype=2)+
    geom_smooth(method = "lm",alpha=0.1, se = T, fill="grey50", colour="black", size=2)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab(variable_x_name)+
    labs(title=paste(variable_x_name,"vs","Richness"))+
    theme_classic()+
    theme(legend.position="none")
}

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=shan,fill=ID),shape=21, size=6)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            theme_classic()+theme(legend.direction = "horizontal",legend.box="vertical"))


png(filename ="Figure/Biol_treat/Shannon.png", 
    width = 1000*4, height = 1000*4, 
    units = "px",res = 300) 
grid.arrange(
  plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],plots_HOB_BDD[[3]],plots_HOB_BDD[[4]], 
  plots_HOB_BDD[[5]],plots_HOB_BDD[[6]],plots_HOB_BDD[[7]],plots_HOB_BDD[[8]], 
  plots_HOB_BDD[[9]],plots_HOB_BDD[[10]],plots_HOB_BDD[[11]],plots_HOB_BDD[[12]], 
  plots_HOB_BDD[[13]],plots_HOB_BDD[[14]],plots_HOB_BDD[[15]],plots_HOB_BDD[[16]], 
  plots_HOB_BDD[[17]],plots_HOB_BDD[[18]],plots_HOB_BDD[[19]], legend_plots,
  nrow=5,ncol=4)
dev.off()

# 5.3. Matching with Bray & Jaccard ####
output <- matrix(nrow=length(unlist(HOB_BDD_match_matrix$BrayCurtis)))
for (lis_elem in 1:length(HOB_BDD_match_matrix)) {
  out <- c()
  for (river in 1:length(HOB_BDD_match_matrix$BrayCurtis)) {
    col_names_river <- names(HOB_BDD_match_matrix$BrayCurtis)
    
    out_temp <- c()
    out_temp <- cbind(as.matrix(HOB_BDD_match_matrix[[lis_elem]][[river]])[
      upper.tri(as.matrix(HOB_BDD_match_matrix[[lis_elem]][[river]]))])
    out_temp <- as.data.frame(cbind("ID"=rep(col_names_river[river],nrow(out_temp)),out_temp))
    out<- rbind(out, out_temp)
  }
  colnames(out)[2] <- names(HOB_BDD_match_matrix)[lis_elem]
  output[,1] <- out[,1]
  output <- as.data.frame(cbind(output, out))
  output <- output%>%select(-ID)
  output[,lis_elem+1] <- as.numeric(output[,lis_elem+1])
}
colnames(output)[1] <- "ID"

STmatrix_BiolDissim <- output

# 5.3.1. Plot Bray ####
plots_HOB_BDD <- list()
for (col_var in 1:ncol(output%>%select(-ID, -BrayCurtis, -Jaccard))){
  
  variable_x <- output%>%select(-ID, -BrayCurtis, -Jaccard)
  variable_x_name <- colnames(variable_x)[col_var]
  variable_x <- variable_x[,col_var]
  
  variable_y <- output%>%select(BrayCurtis)
  
  factors <- output%>%select(ID)
  
  dataset <- cbind(factors, "X_var"=variable_x, variable_y)
  plots_HOB_BDD[[col_var]] <- ggplot(dataset,aes(x=X_var, y=BrayCurtis,fill=ID))+
    geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
    geom_smooth(aes(colour=ID),method = "lm",alpha=0.1, se = F, linetype=2)+
    geom_smooth(method = "lm",alpha=0.1, se = T, fill="grey50", colour="black", size=2)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab(variable_x_name)+
    labs(title=paste(variable_x_name,"vs","Richness"))+
    theme_classic()+
    theme(legend.position="none")
}

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=BrayCurtis,fill=ID),shape=21, size=6)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))

png(filename ="Figure/Biol_treat/Bray.png", 
    width = 1250*2, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],legend_plots,plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
  nrow=2,ncol=3)
dev.off()

# 5.3.2. Plot Jaccard ####
plots_HOB_BDD <- list()
for (col_var in 1:ncol(output%>%select(-ID, -BrayCurtis, -Jaccard))){
  
  variable_x <- output%>%select(-ID, -BrayCurtis, -Jaccard)
  variable_x_name <- colnames(variable_x)[col_var]
  variable_x <- variable_x[,col_var]
  
  variable_y <- output%>%select(Jaccard)
  
  factors <- output%>%select(ID)
  
  dataset <- cbind(factors, "X_var"=variable_x, variable_y)
  plots_HOB_BDD[[col_var]] <- ggplot(dataset,aes(x=X_var, y=Jaccard,fill=ID))+
    geom_point(color="grey30", alpha=0.5, shape=21, size=2)+
    geom_smooth(aes(colour=ID),method = "lm",alpha=0.1, se = F, linetype=2)+
    geom_smooth(method = "lm",alpha=0.1, se = T, fill="grey50", colour="black", size=2)+
    scale_fill_CUNILLERA(palette = "LGTBI")+
    scale_color_CUNILLERA(palette = "LGTBI")+
    xlab(variable_x_name)+
    labs(title=paste(variable_x_name,"vs","Richness"))+
    theme_classic()+
    theme(legend.position="none")
}

legend_plots<- get_legend(ggplot(dataset)+
                            geom_point(aes(x=X_var,y=Jaccard,fill=ID),shape=21, size=6)+
                            scale_fill_CUNILLERA(palette = "LGTBI", name="Stream ID")+
                            theme_classic()+theme(legend.direction = "vertical",legend.box="vertical"))

png(filename ="Figure/Biol_treat/Jaccard.png", 
    width = 1250*2, height = 1000*2, 
    units = "px",res = 300) 
grid.arrange(plots_HOB_BDD[[1]],plots_HOB_BDD[[2]],legend_plots,plots_HOB_BDD[[3]],plots_HOB_BDD[[4]],
             nrow=2,ncol=3)
dev.off()


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# 6. FINAL RESULTS TABLES ####
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#Matrices in list format and in datasframe 
HOB_BDD_match_matrix
STmatrix_BiolDissim
#Individual values altogether 
HOB_BDD_match

