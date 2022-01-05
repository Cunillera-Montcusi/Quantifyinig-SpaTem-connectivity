# Load data management packages
library(dplyr)
library(tidyr)
library(ggplot2)
# Upload packages
library(nlme)
library(lsmeans)
library(multcomp)
library(car)
library(multcompView)

setwd("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/1. FEHM coses al DIA/4. Mecodispers Spa-Tem/MECODISPER SpaTem")

# Data uploading
BDD <- read.csv2("BiolData/Matriz_otoño.csv", sep=";")

BDD[which(BDD$Code=="SA4")[1],3:ncol(BDD)] <- apply(BDD[which(BDD$Code=="SA4"),3:ncol(BDD)],2,sum)
BDD <- BDD[-which(BDD$Code=="SA4")[2],]

Dades_Tardor <- BDD %>%
                select(-Date)%>% # Eliminamos la columna 1 que no tiene nada (select con un "-" = eliminamos)
                pivot_longer(-Code)

colnames(Dades_Tardor) <- c("Code","Species","Abund")

#RICHNESS ________________________________________________________________________####
by_group_dades_tardor <- group_by(Dades_Tardor, 
                                  Code, Species)

data_G_gamma <- by_group_dades_tardor%>%
                group_by(Species)%>%
                summarise(n=sum(Abund))%>%
                mutate(s=1)%>%
                summarise(rich=sum(s))

data_A_alpha <- by_group_dades_tardor%>%
                group_by(Code)%>%
                mutate(s=ifelse(Abund>0,1,0))%>%
                summarise(rich=sum(s))
                
library(vegan)
table_shannon <- by_group_dades_tardor%>%
    spread(key =Species,value = Abund, fill=0)
    
shan <- diversity(table_shannon[,-1],index = "shannon")


Bray.distance<-vegdist(
              decostand(table_shannon[,-1],method = "hellinger"),
              method = "bray") #Càlcul de la matriu de similitud amb Bray-Curtis

table_shannon_PA <- by_group_dades_tardorardor%>%
  mutate(PA=ifelse(Abund>0,1,0))%>%
  select(-Abund)%>%
  spread(key =Species,value = PA, fill=0)

Jac.distance<-vegdist(
              decostand(table_shannon_PA[,-1],method = "pa"),
              method = "jaccard") #Càlcul de la matriu de similitud amb Bray-Curtis


BID_output_IND <- data.frame(Code=unique(by_group_dades_tardor$Code),rich=data_A_alpha$rich,shan)

BID_output_DIST <- list(Code=unique(by_group_dades_tardor$Code),Bray.distance,Jac.distance)

SampSites <- read.table("BiolData/Sites.txt",sep = ",",header = T, dec=".")






