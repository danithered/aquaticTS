##########################################################################################
##########            Seasonality Project, a discrete-time model
##########  Evolutionary Ecology Project
##########  Andras Csercsa, Samuel Dijoux
##########################################################################################
# setwd( "C:/Users/dijous00/Dropbox/Evolutionary-Ecology_Project")
load()
########## Model Parameter =======
## Parameters that will be used for the model simulation
rho = 0.3   #or 0.2       # Flow-through rate #day-1
K   = 3E2   # Resource Carying capacity       #g.L-1
eps = 0.50  # Consumer energy conversion efficiency    #constant
mu_1= 0.09  # Natural mortality of consumer 1 #day-1
mu_2= 0.09  # Natural mortality of consumer 2 #day-1
#Alpha_I = 50  # Species peak performance

########## Temperature gradient setting =======
Tmin_simul <- 0     # Minimal temperature treshold for the simulation
Tmax_simul <- 35    # Maximal temperature treshold for the simulation

D_T = Tmax_simul-Tmin_simul             # Temperature range
n_temp_values = 1000                    # Number of temperature values
## DB: why so many?
## SD: This enable to smooth species temperature-dependence performance distribution (see below) 

T_grad <- seq(Tmin_simul,Tmax_simul,length=n_temp_values) # Temperature gradient

########## Time length of the simulation =======
## Creation of the time sequence

#Times <- seq(0, 363, length = 364)                # Year
Times <- seq(0, 364*3, length = 364*3) ## 3 year simulation

########## Species Thermal niche and temperature-dependence performance =======
N_comp <- 2                              # Number of competitive species

## Creation of 3 thermal niches matrices that record species minimal and maximal suitable 
### Scenario I: Small overlap of Species thermal niches
ThN_I <- matrix(ncol=N_comp, nrow=2)  # Species thermal niche matrix
ThN_I[,1] <- c(8,23)                  # Sp.1 is present from 8-22 degrees C
ThN_I[,2] <- c(15,30)                 # Sp.2 is present from 15-30 degrees C

### Scenario II: High overlap of Species thermal niches
ThN_II <- matrix(ncol=N_comp, nrow=2)  # Species thermal niche matrix
ThN_II[,1] <- c(10,25)                 # Sp.1 is present from 10-25 degrees C
ThN_II[,2] <- c(12,27)                 # Sp.2 is present from 12-27 degrees C

### Scenario III: Generalist vs. Specialist
ThN_III <- matrix(ncol=N_comp, nrow=2)  # Species thermal niche matrix
ThN_III[,1] <- c(15,23)                 # Sp.1 is specialist, present only from 15 to 23 degrees C
ThN_III[,2] <- c(5,30)                  # Sp.2 is generalist, present for 5 to 30 degrees C

## Creation of the species presence/absence matrices along the temperature gradient
### regarding their thermal niche for each scenario
# Null matrices
S_ThN_I   <- matrix(0, nrow=length(T_grad),ncol=N_comp) 
S_ThN_II  <- matrix(0, nrow=length(T_grad),ncol=N_comp) 
S_ThN_III <- matrix(0, nrow=length(T_grad),ncol=N_comp)

# Simulation to replace the null values by species suitable conditions 
for(i in 1:N_comp){
  for(j in 1:length(T_grad)){
    if(T_grad[j]<ThN_I[1,i] || T_grad[j]>ThN_I[2,i]){     # for scenario I
      S_ThN_I[j,i]<-0}else{S_ThN_I[j,i]<-1}
    if(T_grad[j]<ThN_II[1,i] || T_grad[j]>ThN_II[2,i]){   # for scenario II
      S_ThN_II[j,i]<-0}else{S_ThN_II[j,i]<-1}
    if(T_grad[j]<ThN_III[1,i] || T_grad[j]>ThN_III[2,i]){ # for scenario III
      S_ThN_III[j,i]<-0}else{S_ThN_III[j,i]<-1}
  }
}

## Species temperature-dependence performance
# We first create Null Matrices that will contain species performance along the thermal gradient
S_Pf_I <- matrix(0 ,nrow=length(T_grad),ncol=N_comp)
S_Pf_II <- matrix(0 ,nrow=length(T_grad),ncol=N_comp)
S_Pf_III <- matrix(0 ,nrow=length(T_grad),ncol=N_comp)

# For each scenario
# We create a Beta density-distribution that will be used for temperature-dependence performance
# We then use a function that match species presence along the thermal gradient with the Beta density distribution to create species temperature-dependence performance

### Sceranio I
Alpha_I = 50 # Multiplicative constant to calculate species performance for scenario I

x_sI <- cbind(seq(0, 1, length = length(S_ThN_I[S_ThN_I[,1]==1,1])),
        seq(0, 1, length = length(S_ThN_I[S_ThN_I[,2]==1,2])))
dx_sI <- dbeta(x_sI, 5, 2) # using a Beta distribution
# 
for(i in 1:N_comp){
  k <- 1
  for(j in 1:dim(S_ThN_I)[1]){
    if(S_ThN_I[j,i] == 0){
      S_Pf_I[j,i] <- 0
    }else{
      S_Pf_I[j,i]<-as.numeric(dx_sI[k,i]*Alpha_I)
      k <- k+1
    }
  }
}

### Sceranio II
Alpha_II = 50 # Multiplicative constant to calculate species performance for scenario II

x_sII <- cbind(seq(0, 1, length = length(S_ThN_II[S_ThN_II[,1]==1,1])),
              seq(0, 1, length = length(S_ThN_II[S_ThN_II[,2]==1,2])))
dx_sII <- dbeta(x_sII, 5, 2) # using a Beta distribution
#
for(i in 1:N_comp){
  k <- 1
  for(j in 1:dim(S_ThN_II)[1]){
    if(S_ThN_II[j,i] == 0){
      S_Pf_II[j,i] <- 0
    }else{
      S_Pf_II[j,i]<-as.numeric(dx_sII[k,i]*Alpha_I)
      k <- k+1
    }
  }
}

### Sceranio III
Alpha_III_1 = 60 # Multiplicative constant for specialist species, scenario III
Alpha_III_2 = 30  # Multiplicative constant for generalist species, scenario III

x_sIII_1 <- seq(0, 1, length = length(S_ThN_III[S_ThN_III[,1]==1,1]))
x_sIII_2 <- seq(0, 1, length = length(S_ThN_III[S_ThN_III[,2]==1,2]))
dx_sIII_1 <- dbeta(x_sIII_1, 5, 2) # using a Beta distribution
dx_sIII_2 <- dbeta(x_sIII_2, 5, 2) # using a Beta distribution
#
for(i in 1:N_comp){
  k <- 1
  if(i == 1){
      D = dx_sIII_1; ALPHA = Alpha_III_1
    }else{
      D = dx_sIII_2; ALPHA = Alpha_III_2
    }
  
  for(j in 1:dim(S_ThN_III)[1]){
    if(S_ThN_III[j,i] == 0){
      S_Pf_III[j,i] <- 0
    }else{
      S_Pf_III[j,i]<-as.numeric(D[k]*ALPHA)
      k <- k+1
    }
  }
}

########## Plot of species thermal niche per scenarii =======
plot.new()
par(mfrow = c(1, 3))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### Scenario I ========================
par(mfg=c(1,1), mar = c(1, 4, 1, 0.5))
plot(S_Pf_I[,1]~T_grad , type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,max(T_grad)), ylim=c(0, 250), lwd=1, col = "blue", main= "Scenario I")
axis(1, at=c(0, 10, 20, 30), labels=c(0, 10, 20, 30))
axis(2, at=c(0,100,200), labels=c(0, 100, 200))
mtext("Temperature", 1, line=2.0, cex=1.1)
mtext("Feeding rate", 2, line=2.0, cex=1.1)
lines(S_Pf_I[,2]~T_grad, col = "red")

### Scenario II ========================
par(mfg=c(1,2), mar = c(1, 4, 1, 0.5))
plot(S_Pf_II[,1]~T_grad , type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,max(T_grad)), ylim=c(0, 250), lwd=1, col = "blue", main= "Scenario II")
axis(1, at=c(0, 10, 20, 30), labels=c(0, 10, 20, 30))
axis(2, at=c(0,100,200), labels=c(0, 100, 200))
mtext("Temperature", 1, line=2.0, cex=1.1)
mtext("Feeding rate", 2, line=2.0, cex=1.1)
lines(S_Pf_II[,2]~T_grad, col = "red")

### Scenario III ========================
par(mfg=c(1,3), mar = c(1, 4, 1, 0.5))
plot(S_Pf_III[,1]~T_grad , type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,max(T_grad)), ylim=c(0, 250), lwd=1, col = "blue", main= "Scenario III")
axis(1, at=c(0, 10, 20, 30), labels=c(0, 10, 20, 30))
axis(2, at=c(0,100,200), labels=c(0, 100, 200))
mtext("Temperature", 1, line=2.0, cex=1.1)
mtext("Feeding rate", 2, line=2.0, cex=1.1)
lines(S_Pf_III[,2]~T_grad, col = "red")


########## Dataset for Thermal niche
Data_nichab <- cbind(T_grad,Sp_Perf[,c(1,2)]); 
Data_nichab <- as.data.frame(Data_nichab)
colnames(Data_nichab) = c("Tgrad","Sp_perf1", "Sp_perf2")

plot(Data_nichab[,2]~Data_nichab[,1],xlab = "Temperature gradient", ylab="Species performance",
     type = "l", col="blue", lwd=2)
lines(Data_nichab[,3]~Data_nichab[,1], col="red", lwd=2)


########## Temperature variation =======
### This step can be skiped if you load "Simul_T.R"

# We simulate a stochastic environment (i.e. daily fluctuation of temperatures) 
T0 <- 15            # Initial temperature at the beginning of the simulation
Temperature <- NULL
Temperature <- T0   # Temperature record

## Simulation
for(i in 2:length(Times)-1){
  T1 <- rnorm(1, mean=T0, sd=2)                  
  T1 <- min(max(T1, Tmin_simul, na.rm=TRUE), Tmax_simul, na.rm=TRUE)
  Temperature <- append(Temperature,T1)
  T0 <- T1
}
Simul_T<-Temperature
save(Simul_T, file = "Simul_T.R")

#### Simulation with set.seed()
## DB: use set.seed() for reproducible research!
## SD: The problem by using use.seed is that the temperature will only either increase to T_max or decrease to T_min 
for(i in 2:length(Times)){
  set.seed(3)
  T1 <- rnorm(1, mean=T0, sd=2)                  
  T1 <- min(max(T1, Tmin_simul, na.rm=TRUE), Tmax_simul, na.rm=TRUE)
  Temperature <- append(Temperature,T1)
  T0 <- T1
}
plot(Simul_T, type='l', xlab="day", ylab="Temperatures")
########## Daily fluctuation of species feeding rates =======
load("Simul_T.R")

## For each scenario, we record species varying feeding rate per day
### Scenario I   =======
P_1<-NULL; P_2<-NULL;p_1<-NULL;  p_2<-NULL

for(i in 1:length(Simul_T)){
  j   <- NULL
  for(j in 1:length(T_grad)-1){
    if((Simul_T[i] == max(min(Simul_T[i],T_grad[j+1]),T_grad[j-1])) || (Simul_T[i] == min(max(Simul_T[i],T_grad[j+1]),T_grad[j-1]))){
       p_1 <- S_Pf_I[j,1]
       p_2 <- S_Pf_I[j,2]
    }
  }
  P_1 <- append(P_1,p_1)
  P_2 <- append(P_2,p_2)
}
Data_SI <- cbind(P_1, P_2)

Data_SI <- as.data.frame(Data_SI)
colnames(Data_SI) = c("Sp1_Fr_I", "Sp2_Fr_I")
names(Data_SI);attach(Data_SI)

## Observation =======
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

abline(h=ThN_I[1,1],lty=3, col="blue", lwd=2)
abline(h=ThN_I[2,1],lty=4, col="blue", lwd=2)
abline(h=ThN_I[1,2],lty=3, col="red", lwd=2)
abline(h=ThN_I[2,2],lty=4, col="red", lwd=2)

### Sp1 ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(Sp1_Fr_I, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "blue")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[1](T), 2, line=3.5, cex=1.5)
### Sp2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(Sp2_Fr_I, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "red", main ="Scenario I")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[2](T), 2, line=3.5, cex=1.5)


### Scenario II  =======
T_C<-NULL; P_1<-NULL; P_2<-NULL;p_1<-NULL;  p_2<-NULL

for(i in 1:length(Simul_T)){
  j   <- NULL
  for(j in 1:length(T_grad)-1){
    if((Simul_T[i] == max(min(Simul_T[i],T_grad[j+1]),T_grad[j-1])) || (Simul_T[i] == min(max(Simul_T[i],T_grad[j+1]),T_grad[j-1]))){
      p_1 <- S_Pf_II[j,1]
      p_2 <- S_Pf_II[j,2]
    }
  }
  P_1 <- append(P_1,p_1)
  P_2 <- append(P_2,p_2)
}
Data_SII <- cbind(P_1, P_2)

Data_SII <- as.data.frame(Data_SII)
colnames(Data_SII) = c("Sp1_Fr_II", "Sp2_Fr_II")
names(Data_SII);attach(Data_SII)


## Observation =======
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

abline(h=ThN_II[1,1],lty=3, col="blue", lwd=2)
abline(h=ThN_II[2,1],lty=4, col="blue", lwd=2)
abline(h=ThN_II[1,2],lty=3, col="red", lwd=2)
abline(h=ThN_II[2,2],lty=4, col="red", lwd=2)

### Sp1 ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(Sp1_Fr_II, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "blue")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[1](T), 2, line=3.5, cex=1.5)
### Sp2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(Sp2_Fr_II, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "red", main = "Scenario II")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[2](T), 2, line=3.5, cex=1.5)


### Scenario III  =======
T_C<-NULL; P_1<-NULL; P_2<-NULL;p_1<-NULL;  p_2<-NULL

for(i in 1:length(Simul_T)){
  j   <- NULL
  for(j in 1:length(T_grad)-1){
    if((Simul_T[i] == max(min(Simul_T[i],T_grad[j+1]),T_grad[j-1])) || (Simul_T[i] == min(max(Simul_T[i],T_grad[j+1]),T_grad[j-1]))){
      p_1 <- S_Pf_III[j,1]
      p_2 <- S_Pf_III[j,2]
    }
  }
  P_1 <- append(P_1,p_1)
  P_2 <- append(P_2,p_2)
}
Data_SIII <- cbind(P_1, P_2)

Data_SIII <- as.data.frame(Data_SIII)
colnames(Data_SIII) = c("Sp1_Fr_III", "Sp2_Fr_III")
names(Data_SIII);attach(Data_SIII)

## Observation =======
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

abline(h=ThN_III[1,1],lty=3, col="blue", lwd=2)
abline(h=ThN_III[2,1],lty=4, col="blue", lwd=2)
abline(h=ThN_III[1,2],lty=3, col="red", lwd=2)
abline(h=ThN_III[2,2],lty=4, col="red", lwd=2)

### Sp1 ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(Sp1_Fr_III, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "blue")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[1](T), 2, line=3.5, cex=1.5)
### Sp2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(Sp2_Fr_III, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 150), lwd=2, col = "red", main="Scenario III")
axis(2, at=c(50, 100), labels=c(50, 100))
mtext(~f[2](T), 2, line=3.5, cex=1.5)


########## 
## Initialisation phase.  =======
y0 <- pop_initial <- c(100, 0, 0)
R_t <- y0[1];
C1_t <- y0[2];
C2_t <- y0[3]  #Initial population densisies 

popdata<-matrix(0, nrow=length(Times),ncol=4)
popdata[1,1:4] <- c(1, y0[1], y0[2], y0[3])

## Simulation 
for(i in 2:length(Times)){
  
  if(R_t<=0){
    sp1_impact <- 0
    sp2_impact <- 0
  }else{
    sp1_impact <- Data_seasonality[i,3]
    sp2_impact <- Data_seasonality[i,4]
  }
  R_T  <- R_t+rho*(K-R_t) - (sp1_impact+sp2_impact) #Pop density at time t=t+1
  C1_T <- C1_t + eps*sp1_impact - mu_1*C1_t
  C2_T <- C2_t + eps*sp2_impact - mu_2*C2_t
  
  R_t <- R_T;
  C1_t <- C1_T;
  C2_t <- C2_T  #Pop density at time t
  popdata[i,1:4] <- c(i,R_T,C1_T,C2_T)
  #popdata[i,1:2] <- c(i,R_T)#,C1_T)#,C2_T)
}
## Scenario I: =======
for(i in 2:length(Times)){
  if(R_t<=0){
    sp1_impact <- 0
    sp2_impact <- 0
  }else{
    sp1_impact <- Data_SI[i,3]
    sp2_impact <- Data_SI[i,4]
  }
  R_T  <- R_t+rho*(K-R_t) - (sp1_impact+sp2_impact) #Pop density at time t=t+1
  C1_T <- C1_t + eps*sp1_impact - mu_1*C1_t
  C2_T <- C2_t + eps*sp2_impact - mu_2*C2_t
  
  R_t <- R_T;
  C1_t <- C1_T;
  C2_t <- C2_T  #Pop density at time t
  popdata[i,1:4] <- c(i,R_T,C1_T,C2_T)
}
popDyn_I <- popdata

## Scenario II: =======
for(i in 2:length(Times)){
  if(R_t<=0){
    sp1_impact <- 0
    sp2_impact <- 0
  }else{
    sp1_impact <- Data_SII[i,3]
    sp2_impact <- Data_SII[i,4]
  }
  R_T  <- R_t+rho*(K-R_t) - (sp1_impact+sp2_impact) #Pop density at time t=t+1
  C1_T <- C1_t + eps*sp1_impact - mu_1*C1_t
  C2_T <- C2_t + eps*sp2_impact - mu_2*C2_t
  
  R_t <- R_T;
  C1_t <- C1_T;
  C2_t <- C2_T  #Pop density at time t
  popdata[i,1:4] <- c(i,R_T,C1_T,C2_T)
}
popDyn_II <- popdata

## Scenario III: =======
for(i in 2:length(Times)){
  if(R_t<=0){
    sp1_impact <- 0
    sp2_impact <- 0
  }else{
    sp1_impact <- Data_SIII[i,3]
    sp2_impact <- Data_SIII[i,4]
  }
  R_T  <- R_t+rho*(K-R_t) - (sp1_impact+sp2_impact) #Pop density at time t=t+1
  C1_T <- C1_t + eps*sp1_impact - mu_1*C1_t
  C2_T <- C2_t + eps*sp2_impact - mu_2*C2_t
  
  R_t <- R_T;
  C1_t <- C1_T;
  C2_t <- C2_T  #Pop density at time t
  popdata[i,1:4] <- c(i,R_T,C1_T,C2_T)
}
popDyn_III <- popdata

########### Plot Scenario I========================
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

########### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

########### R ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(popDyn_I[,2]~popDyn_I[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 350), lwd=2, col = "darkgreen")
axis(2, at=c(150, 300), labels=c(150, 300))
mtext("Resources", 2, line=3.5, cex=1.5)

########### C1 & C2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(popDyn_I[,3]~popDyn_I[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n", 
     xlim=c(0,length(Times)), ylim=c(0, 450), lwd=2, col="blue", main= "Scenario I")
axis(2, at=c(200, 400), labels=c(200, 400))
mtext("Consumers", 2, line=3.5, cex=1.5)
lines(popDyn_I[,4], type="l", lwd=2, col="red")

########### Plot Scenario II========================
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

########### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

########### R ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(popDyn_II[,2]~popDyn_II[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 350), lwd=2, col = "darkgreen")
axis(2, at=c(150, 300), labels=c(150, 300))
mtext("Resources", 2, line=3.5, cex=1.5)

########### C1 & C2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(popDyn_II[,3]~popDyn_II[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n", 
     xlim=c(0,length(Times)), ylim=c(0, 450), lwd=2, col="blue", main= "Scenario II")
axis(2, at=c(200, 400), labels=c(200, 400))
mtext("Consumers", 2, line=3.5, cex=1.5)
lines(popDyn_II[,4], type="l", lwd=2, col="red")

########### Plot Scenario III========================
plot.new()
par(mfrow = c(3, 1))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

########### T ========================
par(mfg=c(3,1), mar = c(1, 4, 0, 0.5))
plot(Simul_T~Times, type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 35), lwd=1, col = "black")
axis(1, at=c(0, 200, 400, 600, 800, 1000), labels=c(0, 200, 400, 600, 800, 1000))
axis(2, at=c(0,15,30), labels=c(0, 15, 30))
mtext("Time", 1, line=1.5, cex=1.5)
mtext("Temperature", 2, line=3.5, cex=1.5)

########### R ========================
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))
plot(popDyn_III[,2]~popDyn_III[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n",
     xlim=c(0,length(Times)), ylim=c(0, 350), lwd=2, col = "darkgreen")
axis(2, at=c(150, 300), labels=c(150, 300))
mtext("Resources", 2, line=3.5, cex=1.5)

########### C1 & C2 ========================
par(mfg=c(1,1), mar = c(0, 4, 1, 0.5))
plot(popDyn_III[,3]~popDyn_III[,1], type="l", xlab="", ylab="", xaxt="n", xaxs="i",yaxs="i", yaxt="n", 
     xlim=c(0,length(Times)), ylim=c(0, 450), lwd=2, col="blue", main= "Scenario III")
axis(2, at=c(200, 400), labels=c(200, 400))
mtext("Consumers", 2, line=3.5, cex=1.5)
lines(popDyn_III[,4], type="l", lwd=2, col="red")
