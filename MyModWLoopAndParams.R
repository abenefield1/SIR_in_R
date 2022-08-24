require(deSolve)
require(ggplot2)
require(reshape2)
library(dplyr)
require(gridExtra)

rm(list=ls())

# Added a new function to create the data for k classes, so here it just creates the parameters
parameters=data.frame(rec=0.2,trans=0.005, death=.1,birth=.1, mut=0.05)

PopulationSize<-1000
#R0<- (parameters$birth * PopulationSize * parameters$trans) / ((parameters$death*PopulationSize) * (parameters$death + parameters$rec))
#R0


# Time points
#time=seq(from=1,to=150,by=1/365)
time=seq(from=1,to=150)

##### detection function ######
# Function to get detection rate with is an exponential decay function scaled by "a"
det = function(k, a = 1, dMax=0){
  return(dMax*exp(-k*a))
}

#### Build data for n_k infection classes ######
# Builds the initial state data vector. Takes total number of infection classes (n_k, which includes I0, e.g. 20 infection classes run from 0 to 19) and then total pop size. This will be distributed such that only S and I0 have individuals, otherwise all 0s.
build_data <- function(n_k, total_pop_size){
  
  state <- numeric(n_k+2)
  names(state)[c(1,n_k+2)] <- c('S','R')
  
  for(i in 1:n_k){
    names(state)[i+1] <- paste0("I",i-1)
  }
  
  state[1] <- total_pop_size-1
  state[2] <- 1
  return(state)
}

initial_states <- build_data(n_k=3,total_pop_size = PopulationSize)
initial_states

Q<-1000 # REMOVE ME - temp pop size

######### SIR Model for k classes #######
sir_modelAW <- function(time, state, parameters){
  
  # Get total pop size by summing over state vector
  N = sum(state)  
  
  # Total number of states including infection, S, and R
  all_states <- length(state)
  
  # Number of potential mutation/infection states (all states minus susceptible and recovered)
  n = all_states - 2
  
  # vector to hold di/dt 
  i_dot = rep(0, n)
  
  # I_n is the number of infected in the last infection/mutation class
  I_n = state[length(state)-1]
  
  # number of recovered individuals
  R = state[length(state)]
  
  # Vector to store deltas for all states, from S through all I's and N
  deltas <- numeric(all_states)
  
  
  #################################################################################
  ### QUESTION:
  #################################################################################
  i = state[2:n] ############## SHOULD THIS BE 2:n-1 b/c I_n is coded separately???
  #################################################################################
  #################################################################################
  
  
  ## S #########################################
  #First Delta is for susceptibles - don't think this needs to be iterated over all classes b/c of the sum_of_rel function? Could be wrong.
  deltas[1] = parameters$birth*Q - (state['S']*(parameters$death + parameters$trans*sum(state[which(names(state)!='S' & names(state)!='R')])))
  
  
  ### I0 ########################################
  # Second delta is I0 which is unique b/c nothing mutates into it. So it gets its own function here, gaining new individuals from S in the first half, then losing individuals through recovery, mutation to I1, death, and detection.
  deltas[2] = parameters$trans*state['S']*state['I0'] - (parameters$rec + parameters$mut + parameters$death + det(k=0))*state['I0']
  
  ### In ########################################
  #The second to last delta is the last infection class, which is also unique, in that it doesn't mutate into anything. So it gets its own recursion. Gaining individuals from S in the first half, losing individuals from recovery, death, and detection, then gaining individuals from the second to last infection class via mutation.
 deltas[all_states-1] = parameters$trans*state['S']*I_n - (parameters$rec + parameters$death + det(k=n))*I_n + parameters$mut*state[all_states-2]
  
  # Pulling the indices for all infection classes except for I0 and In, since they should all have the same equation for change over time: gains from susceptibles based on transmission, loss from recovery, mutation, death, and detection, then gain from mutation in the previous class.
  
  # Construct the last infection class name using n
  last_name <- paste0('I', n)
  # vector of names to exclude
  excluded_names <- c("S", "I0", last_name, "R")
  
  #indices for all infection classes. In retrospect, this should always be like 3:all_states-2 I think? We exclude S and I0 (indices 1 and 2) and then In and R (indices all_states-1 and all_states)
  inf_indices <- which((!names(state) %in% excluded_names))
  
  # Not sure If I need to do this, but I've set this up do recovery from class I0 and class In as well as death of recovered individuals. I then update this with recoveries from all other infection classes within the loop below. Then I put that into the deltas vector at the end. there may be a more elegant way of doing this using the sum_of_rel variable/function but I wasn't sure how that worked.
  dR <- (parameters$rec +det(k=0))*state['I0']+(parameters$rec + det(k=all_states-2))*state[all_states-1] - parameters$death*state['R']
  
  # This iterator tracks how many mutations we are from 0, which is just tough to do when using indices, but juggling indices is tough if we try and iterate over the infection class number itself. So I just got lazy about it tbh. Hopefully it works?
  it <- 1
  
  for(i in inf_indices){
    
    # Get delta for current infection class, looping through all except 0 and n
    deltas[i] <- parameters$trans*state['S']*state[i] - (parameters$rec + parameters$mut + parameters$death + det(k=it))*state[i] + parameters$mut*state[i-1]
    
    #get our recovered change - just add to it each time I think? This makes sense? Maybe? Hmmm.
    dR <- dR + (parameters$rec + det(k=it))*state[i] - parameters$death*state['R']
    
    #update our mutation iterator
    it <- it+1
  }
  
  #dR <- dR - parameters$death*(state['R']+dR)
  
  #   #put the recovered change into the final deltas slot
  deltas[all_states] <- dR
  #   
  #   #return deltas
  return(list(deltas))
}

initial_states

sir_modelAW(time, initial_states, parameters)


###### Solving the differential equations: ######
output<-as.data.frame(ode(y=initial_states,func = sir_modelAW,parms=parameters,times = time, method="ode45"))
max(rowSums(output[,-1]))

# Get the sums of all infected classes to plot as a separate line
sums<-output[,c(-2, -length(output))] # remove the S and R class

Totals<-rowSums(sums[,-1]) # Take the row sums excluding the time variable

# Format it to fit with out_long
Tots<-data.frame("time"=as.numeric(sums$time), "variable"=rep("Total_I", length(Totals)), "value"=as.numeric(Totals))

out_long=melt(output,id="time")

which(out_long$variable!="S" & out_long$variable!= "I0"  & out_long$variable!= "R")

out_long<-rbind(out_long, Tots)

### Creating a class to relabel for plotting
out_long$class<-NA

out_long$class[which(out_long$variable!="S" & out_long$variable!= "I0"  & out_long$variable!= "R" & out_long$variable!= "Total_I")]<-"I_{k}"
out_long$class[which(out_long$variable=="S")]<-"S"
out_long$class[which(out_long$variable=="I0")]<-"I_{0}"
out_long$class[which(out_long$variable=="R")]<-"R"
out_long$class[which(out_long$variable=="Total_I")]<-"Total_I"

head(out_long)

####### Plotting the prevalence over time #######
#quartz(width=10,height=6)
ggplot(data = out_long,                                               
       aes(x = time, y = value, colour = class, group = variable, size=class)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                  
  ylab("Prevalence") + scale_color_discrete(name="State")+
  scale_color_manual(values=c("red", "#e4ab18", "blue", "darkgreen", "black"))+
  scale_size_manual(values=c(0.5,0.5,0.75,0.75,0.75))

#rm(list=ls())
