# Script for generating and saving scrambled hydrographs, along with their KGE and NSE scores
# This version also tries adding the ensemble for better efficiency
# The bias shift is the same as v1 here

#Load the 'scramble_simgen' function FIRST!

library(hydroGOF)

# This version (Oct22, 2024) the function does the power stretch, drift and smooth, but not shift or multiply
#   (to avoid loading big sim files)
# Need to: run the 'obs' hyd through the function, then nested loops to multiply and add the shift, 
#  Finally calc the stats
# overwrite the 'sim' to keep space use under control

#Obs should be a vector of 365 value (1-year daily flows)
# REPLACE TO RUN
obs<-Var_b42

#Run the scrambling function
sim_base<-scramble_simgen(obs)
sim<-sim_base
n<-ncol(sim_base)
id<-c(1:ncol(sim))
obs_mean<-mean(obs)

#init for ensemble stuff
#Set the thresholds to check (fewer is faster but reduces resolution)
threshold<-c(seq(0,1,by=0.02))
min_NSE<-matrix(1000,365,length(threshold))
max_NSE<-matrix(-1,365,length(threshold))
avg_NSE<-matrix(0,365,length(threshold))
min_KGE<-matrix(1000,365,length(threshold))
max_KGE<-matrix(-1,365,length(threshold))
avg_KGE<-matrix(0,365,length(threshold))
NSE_bandwidth<-matrix(0,length(threshold))
KGE_bandwidth<-matrix(0,length(threshold))
NSE_total<-matrix(0,length(threshold))
KGE_total<-matrix(0,length(threshold))

#Now the nested loops for multiplying then inside, adding
# Adding +/-% of the mean flow: 2,5,10,15,20,25,30,40,50,60,80,100,120,140,160,200,250,300,350,400
add_vector<-c(0,2,5,10,15,20,25,30,40,50,60,80,100,120,140,160,200,250,300,350,400,-2,-5,-10,-15,-20,-25,-30,-40,-50,-60,-80,-100,-120,-140,-160,-200,-250,-300,-350,-400)/100
# Multiply by 1+/-: 0.01,0.02,0.04,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,0.999
mult_vector<-1+c(0,0.01,0.02,0.04,0.06,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,0.999,-0.01,-0.02,-0.04,-0.06,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5,-0.6,-0.7,-0.8,-0.9,-0.999)

#Outer loop, multiply:
for(m in 1:length(mult_vector)){
  print(m)
  #Inner loop: add
  for(p in 1:length(add_vector)){
    print(p)
    #Adjust the simulation
    sim<-mult_vector[m]*sim_base+add_vector[p]*obs_mean
    sim<-replace(sim,sim<0,0)
    hyd_temp<-sim
    #initilize the stat storage matrixes
    KGE_set<-matrix(0,n)
    NSE_set<-matrix(0,n)
    #calculate the stats for every simulation
    for(i in 1:n){
      KGE_set[i]<-KGE(sim[1:365,i],obs)
      NSE_set[i]<-NSE(sim[1:365,i],obs)
    }
 
        #Now the ensemble stuff:
    NSE_set<-replace(NSE_set,is.na(NSE_set),-99)
    KGE_set<-replace(KGE_set,is.na(KGE_set),-99)
    hyd_df<-data.frame(KGE=KGE_set,NSE=NSE_set,ID=id)
    #With the small section ready, check if it changes the max or min at all threshold levels
    for(j in 1:length(threshold)){
      #count all the sims in the section over the threshold
      NSE_count<-length(NSE_set[NSE_set>=threshold[j]])
      KGE_count<-length(KGE_set[KGE_set>=threshold[j]])
      #pull all the values in the section over the threshold (to get the ids)
      NSE_df<-subset(hyd_df,NSE>=threshold[j])
      KGE_df<-subset(hyd_df,KGE>=threshold[j])
      #If there are any sims over the threshold...
      if(NSE_count>0){
        #Make a sub-set matrix
        subset_sim<-matrix(hyd_temp[1:365,NSE_df$ID],365,NSE_count)
        NSE_total[j]<-NSE_total[j]+NSE_count
        #Then check if it changes the outer limits for a threshold and update if so
        for(t in 1:365){
          if(min(subset_sim[t,1:NSE_count])<min_NSE[t,j]){min_NSE[t,j]<-min(subset_sim[t,1:NSE_count])}
          if(max(subset_sim[t,1:NSE_count])>max_NSE[t,j]){max_NSE[t,j]<-max(subset_sim[t,1:NSE_count])}
          avg_NSE[t,j]<-avg_NSE[t,j]+sum(subset_sim[t,1:NSE_count])
        }}
      
      #If there are any sims over the threshold...
      if(KGE_count>0){
        #Make a sub-set matrix
        subset_sim<-matrix(hyd_temp[1:365,KGE_df$ID],365,KGE_count)
        KGE_total[j]<-KGE_total[j]+KGE_count
        #Then check if it changes the outer limits for a threshold and update if so
        for(t in 1:365){
          if(min(subset_sim[t,1:KGE_count])<min_KGE[t,j]){min_KGE[t,j]<-min(subset_sim[t,1:KGE_count])}
          if(max(subset_sim[t,1:KGE_count])>max_KGE[t,j]){max_KGE[t,j]<-max(subset_sim[t,1:KGE_count])}
          avg_KGE[t,j]<-avg_KGE[t,j]+sum(subset_sim[t,1:KGE_count])
        }}
    }
  }
}

#After the global max/min are found for all thresholds, calculate:
for(j in 1:length(threshold)){  
  KGE_bandwidth[j]<-mean((max_KGE[1:365,j]-min_KGE[1:365,j])/obs)
  NSE_bandwidth[j]<-mean((max_NSE[1:365,j]-min_NSE[1:365,j])/obs)
  avg_KGE[1:365,j]<-avg_KGE[1:365,j]/KGE_total[j]
  avg_NSE[1:365,j]<-avg_NSE[1:365,j]/NSE_total[j]
}

b42_bandwidth<-data.frame(Threshold=threshold,KGE=KGE_bandwidth,NSE=NSE_bandwidth)
b42_NSEmax<-max_NSE
b42_NSEmin<-min_NSE
b42_KGEmax<-max_KGE
b42_KGEmin<-min_KGE
b42_KGEavg<-avg_KGE
b42_NSEavg<-avg_NSE