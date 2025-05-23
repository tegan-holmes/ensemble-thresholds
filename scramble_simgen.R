scramble_simgen<- function(obs){
  # This function takes an 'observed' hydrograph and scramble it progressively
  #to generate increasingly wrong 'simulated' hydrographs
  # The output is a massive block of time series, intended to test threshold choices
  # The 'obs' should be 1 year of daily flows
  # **Code is definitely not optimized
  
  #General procedure is to modify, then append, such that the 'sim' multiply
  #sim<-cbind(sim,obs)
  
  #This version scrambles in 3 error dimensions: smoothing, powers and temporal shifting
  #Added and multiplied bias is done later to the saved scrambled set
  
  #init:
  sim<-obs
  
  #Smoothing: running averages
  # The intervals:
  smooth_increments<-c(2,3,5,7,10,15,20,25,30,45,60,90,120,150,180,210,300)
  
  for(j in 1:length(smooth_increments)){
    n<-smooth_increments[j]
    #initilize
    sim_smooth<-obs
    #loop through each day
    for(t in 1:365){
      start_day<-t-floor((n-1)/2)
      end_day<-t+ceiling((n-1)/2)
      if(start_day<1){Q_sum<-sum(obs[(365+start_day):365])+sum(obs[1:end_day])}
      else if(end_day>365){Q_sum<-sum(obs[start_day:365])+sum(obs[1:(end_day-365)])}
      else {Q_sum<-sum(obs[start_day:end_day])}
      sim_smooth[t]<-Q_sum/n
    }
    sim<-cbind(sim,sim_smooth)
  }

  
  #Stretching: multiply by 0.5 to 1.5, in 0.05 increments (skip x1 as it is already there)
#  sim_base<-sim #set the current base case
#  for (stretch_factor in c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5)){
#    sim_temp<-sim_base*stretch_factor
#    sim<-cbind(sim,sim_temp)
#  }
  
  #Stretching: exponents 0.75 to 1.2, in 0.05 increments (skip x1 as it is already there as base case)
  sim_base<-sim #set the current base case [obs+smoothed]
  stretch_factor<-c(0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.02,1.04,1.06,1.08,1.1,1.15,1.2,1.25)
  for(j in 1:length(stretch_factor)){
    sim_temp<-sim_base^stretch_factor[j]
    sim<-cbind(sim,sim_temp)
  }
  
#  #Shifting: add -5 to 5, in 1 increments (skip x1 as it is already there)
#  sim_base<-sim #set the current base case
#  for (shift_factor in c(-5,-4,-3,-2,-1,1,2,3,4,5)){
#    sim_temp<-sim_base+shift_factor
#    sim_temp<-replace(sim_temp,sim_temp<0,0)
#    sim<-cbind(sim,sim_temp)
#  }
  
  #Drifting: move the whole year, in 1 increments (skip x1 as it is already there)
  sim_base<-sim #set the current base case [(obs+smoothed)*powers]
  sim_length<-ncol(sim) #current number of sims
  for (drift_num in c(1:364)){
    sim_temp<-sim_base
    sim_temp[1:drift_num,1:sim_length]<-sim_base[(366-drift_num):365,1:sim_length]
    sim_temp[(drift_num+1):365,1:sim_length]<-sim_base[1:(365-drift_num),1:sim_length]
    sim<-cbind(sim,sim_temp)
  }
  
  return(sim)}
