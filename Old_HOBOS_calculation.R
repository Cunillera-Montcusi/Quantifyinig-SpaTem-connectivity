#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Previous HOBOS values__________________________________________________________________ ####

# Needed values from the script "SpaTemp_HOBOS_treatment.R"
## You need to run the previous script in order to built the HOBOS dataset (a list with the HOBOS info for each river)
HOBOS_sites
# This two other objects are also built in the "SpaTemp_HOBOS_treatment.R" script. They contain the names 
## of each river and the order upstream to downstream within it. 
HOB_riv_ID
ups_dos

## CALCULATION OF THE HOBOS VALUES (as before)
# The total duration of drying events (TotDur) which represents the total number of zero-flow days
# The frequency of drying events (TotNum), which represents the number of drying events
# Mean length (in days) of drying events (TotLeng). 

Old_HOBOS <- list()
TotDur <- list()
TotNum <- list()
TotMean_duration <- list()
for (river in 1:length(HOBOS_sites)) {
  out <- c()
  out_TotNum <- c()
  out_TotNum_rew <- c()
  #TotDur calculation
  cols_repeat <- ncol(HOBOS_sites[[river]])-1
  for (colu in 1:cols_repeat) {
    out[colu] <- length(which(HOBOS_sites[[river]][,colu+1]==0))
    
    #TotNum calculation
    out_sum <- c()
    for (place in 1:(length(HOBOS_sites[[river]][,colu+1])-1)) {
      out_sum[place] <- HOBOS_sites[[river]][,colu+1][place]-HOBOS_sites[[river]][,colu+1][place+1]
    }
    out_TotNum[colu] <- length(which(out_sum==1))
    #Conditions to count the TotMean_duration - The mean duration of each dry period. 
    if(length(which(out_sum==1))-length(which(out_sum==-1))<0){
      out_TotNum_rew[colu] <- mean(which(c(1,out_sum)==1)-which(c(1,out_sum)==-1))}
    if(length(which(out_sum==1))-length(which(out_sum==-1))>0){
      out_TotNum_rew[colu] <- mean(which(c(out_sum,-1)==1)-which(c(out_sum,-1)==-1))}
    if(length(which(out_sum==1))-length(which(out_sum==-1))==0){
      if(out_sum[1]==-1){
        if(length(which(c(1,out_sum)==1))-length(which(c(1,out_sum)==-1))>0){
          out_TotNum_rew[colu] <- mean(which(c(1,out_sum,-1)==1)-which(c(1,out_sum,-1)==-1))
        }}else{
          out_TotNum_rew[colu] <- mean(which(out_sum==1)-which(out_sum==-1))}
    }
  }
  TotDur[[river]] <- out
  TotNum[[river]] <- out_TotNum
  TotMean_duration[[river]] <- -out_TotNum_rew
  TotMean_duration[[river]][is.na(TotMean_duration[[river]])] <- 0
  
  Old_HOBOS[[river]]<- data.frame("TotDur"=TotDur[[river]],
                                  "TotNum"=TotNum[[river]],
                                  "TotLeng"=TotMean_duration[[river]])
}

# Adding ID and DtoU in Old HOBOS values to then match with ST indices
Old_HOBOS
Old_HOBOS_comb <- data.frame("ID"=HOB_riv_ID,"DtoU"=ups_dos,
                             rbind(Old_HOBOS[[1]],Old_HOBOS[[2]],
                                   Old_HOBOS[[3]],Old_HOBOS[[4]],
                                   Old_HOBOS[[5]],Old_HOBOS[[6]],
                                   Old_HOBOS[[7]]))%>%
                             mutate(ID_UpDo=paste(ID,"_",DtoU,sep = ""))