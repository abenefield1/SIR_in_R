require(deSolve)
require(ggplot2)
require(reshape2)
library(dplyr)
require(gridExtra)


# Added a new function to create the data for k classes, so here it require(deSolve)
require(ggplot2)
require(reshape2)
require(gridExtra)

# Function to get detection rate with is an exponential decay function scaled by "a"
det = function(k, a = 0.1, dMax=0){
  return(dMax*exp(-k*a))
}
# Test det(k=2,a=0.5)

initial_state_values=c(S=999,I0=1, I1=0, I2=0,R=0)

# Parameters
parameters=c(rec=0.2,trans=0.005, death=0.1,birth=0.1, mut=0.05, dMax=0, a=0.9)

# Time points

time=seq(from=1,to=150)

# SIR model function

sir_model2 <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S + I0 + I1 + I2 + R
    dS = birth*1000 - (S*(death + trans*I0 + trans*I1 +trans* I2))
    dI0 = trans*S*I0 - (rec + mut + death + dMax*exp(-0*a))*I0
    dI1 = trans*S*I1 - (rec + mut + death + dMax*exp(-1*a))*I1 + mut*I0
    dI2 = trans*S*I2 - (rec + death + dMax*exp(-2*a))*I2 + mut*I1
    dR = ((rec + dMax*exp(-0*a))*I0 + (rec + dMax*exp(-1*a))*I1 + (rec + dMax*exp(-2*a))*I2) - death*R
    
    return(list(c(dS,dI0, dI1, dI2,dR)))
  }
  )
}


# Solving the differential equations:
output<-as.data.frame(ode(y=initial_state_values,func = sir_model2,parms=parameters,times = time, method="ode45"))


# Summing total infected classes:
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

#out_long$class[which(out_long$variable!="S" & out_long$variable!= "I0"  & out_long$variable!= "R" & out_long$variable!= "Total_I")]<-"I_{k}"
out_long$class[which(out_long$variable=="S")]<-"S"
out_long$class[which(out_long$variable=="I0")]<-"I0"
out_long$class[which(out_long$variable=="I1")]<-"I1"
out_long$class[which(out_long$variable=="I2")]<-"I2"
out_long$class[which(out_long$variable=="R")]<-"R"
out_long$class[which(out_long$variable=="Total_I")]<-"Total I"

head(out_long)

#################################################


#Plotting the prevelance over time
ggplot(data = out_long,                                               
       aes(x = time, y = value, colour = class, group = class)) +  
  geom_line() +                                                          
  xlab("Time (years)")+                                                  
  ylab("Prevalence") + scale_color_discrete(name="State") +
  scale_color_manual(values=c("red", "#e4ab18", "#e4ab18", "blue", "darkgreen", "black"))+
  scale_size_manual(values=c(0.5,0.5, 0.5,0.75,0.75,0.75))

  

#rm(list=ls())