#### SEI epidemic model for Celia's water stress study
#### 2015-09-10

rm(list = ls())
setwd("C:/Users/Adam/Documents/UC Berkeley post doc/Almeida lab/Celias water stress study")
my.packages <- c("deSolve", "lattice", "tidyr", "dplyr", "data.table")
lapply(my.packages, require, character.only = TRUE)

#### Functions to produce error bars in lattice plot
source("C:/Users/Adam/Documents/R/Functions/Error_bars_in_lattice_plots.r")


##########################################################################################
#### Numerical simulation
##########################################################################################
# Equations for Water Stress SECI model
# Exactly the same equations as Matt's climate spread model
WSModel <- function(Time, State, Pars){
  with(as.list(c(State, Pars)), {
    dS <- mu*E + mu*C + mu*I - (beta*S*V)/(p*I + S + E + C)
    dE <- (beta*S*V)/(p*I + S + E + C) - gamma*E - mu*E
    dC <- gamma*E - delta*C - mu*C
    dI <- delta*C - mu*I
    dU <- nu*V - (alpha*C*U)/(p*I + S + E + C) - (alpha*p*I*U)/(p*I + S + E + C)
    dV <- (alpha*C*U)/(p*I + S + E + C) + (alpha*p*I*U)/(p*I + S + E + C) - nu*V
    return(list(c(dS, dE, dC, dI, dU, dV)))
  })
}


#### Function to return all the model dynamics
WSdynamicsFunc <- function(x){
  x <- as.numeric(x)
  Pars <- c(alpha = x[1], beta = x[2], delta = x[3], 
            gamma = x[4], p = x[5], mu = x[6], nu = x[7])
  State <- c(S = S, E = E, C = C, I = I, U = U, V = V)
  Time <- seq(0, 2000, by = 1) # time steps to 2000
  model.out <- as.data.frame(ode(func = WSModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  return(model.out)
}


## Function to run simulations and output final (equilibrial) infective host (I) and infectious vector (V) densities
WSsimFunc1 <- function(x){
  x <- as.numeric(x)
  Pars <- c(alpha = x[1], beta = x[2], delta = x[3], 
            gamma = x[4], p = x[5], mu = x[6], nu = x[7])
  State <- c(S = S, E = E, C = C, I = I, U = U, V = V)
  Time <- seq(0, 2000, by = 1) # time steps to 2000 by 1
  model.out <- as.data.frame(ode(func = WSModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  # Need to remove some time steps because under/over flow results in NAs
  model.nona <- model.out[!is.na(model.out$I),]
  model.dat <- data.frame(model.nona[nrow(model.nona),c("time", "C", "I", "U", "V")])
  return(model.dat)
}


## Function to run simulations and output time to equilibrium
WSTimeEqFunc <- function(x){
  x <- as.numeric(x)
  Pars <- c(alpha = x[1], beta = x[2], delta = x[3], 
            gamma = x[4], p = x[5], mu = x[6], nu = x[7])
  State <- c(S = S, E = E, C = C, I = I, U = U, V = V)
  Time <- seq(0, 2000, by = 1) # time steps to 2000
  model.out <- as.data.frame(ode(func = WSModel,
                                 y = State,
                                 parms = Pars,
                                 times = Time))
  # Need to remove some time steps because under/over flow results in NAs
  model.nona <- model.out[!is.na(model.out$I),]
  # Extract time to equilibrium for I
  maxI <- model.nona[nrow(model.nona),"I"] # I* at last time step
  timeEqI <- min(model.nona$time[round(model.nona$I, digits = 5) == round(maxI, digits = 5)])
  # Extract time to equilibtrium for V
  maxV <- model.nona[nrow(model.nona),"V"] # V* at last time step
  timeEqV <- min(model.nona$time[round(model.nona$V, digits = 5) == round(maxV, digits = 5)])
  return(c(maxI, timeEqI, maxV, timeEqV))
}

#### R0 function. Equation from Matt's climate spread ms
R0func <- function(x){
  alpha <- x[1]; beta <- x[2]; delta <- x[3]
  gamma <- x[4]; p <- x[5]; mu <- x[6]; nu <- x[7]
  N <- N; M <- M # Make sure N and M are specified in the R environment
  R0 <- ((alpha*beta*gamma*M*(mu + delta*p))/sqrt(alpha*beta*gamma*mu*nu*N*M*(gamma + mu)*(delta + mu)*(mu + delta*p)))
  return(R0)
}

###################################################################
#### Generate parameter values from experimental data
# Generate nsim simulations, reduce number to nmin, 

# Initial state variables
S <- 100; I <- 0; E <- 0; C <- 0
N <- S + I + E + C
U <- 199; V <- 1
M <- U + V


#### Function to generate simulated data and remove negative (nonsense) values
# Assume normal distribution
simulateData <- function(mean, SE, nsim = 10000, nmin = 1000){
  # nsim = number of simulations
  # nmin = minimum number of final simulations,
  # nmin used to eliminate negative values and maintain symmetry of normal distributions
  rawSim <- rnorm(nsim, mean = mean, sd = SE) # vector of parameter estimates
  simPos <- rawSim[rawSim >= 0] # remove estimates < 0
  simDiff <- length(rawSim) - length(simPos)
  # Remove an equal number of values in the right tail of the distribution
  simVec <- simPos[-order(simPos)[(length(simPos)-simDiff):length(simPos)]][1:nmin]
  return(simVec)
}


#### Unvarying parameters
mu <- 0.05 # Recovery rate of hosts
nu <- 1.5 # Recovery/turnover rate of vectors
gamma <- 0.5
p <- 6/(24-6) # Preference; calculated from Daugherty et al. (2011)
pse <- sqrt((p*(1 - p))/24)  # preference; probability; mean & SE -- MC normal distr
# Generate simulated data
pVec <- simulateData(mean = p, SE = pse)

#### Acquisition and inoculation estimates for Low water stress
# Acquisition
alphalm <- 0.491; alphalse <- 0.066 # acquisition; probability; mean & SE -- MC normal distr
# Generate simulated data
alphalVec <- simulateData(mean = alphalm, SE = alphalse)

# Inoculation
betalm <- 0.360; betalse <- 0.068  # inoculation; probability; mean & SE -- MC normal distr
# Generate simulated data
betalVec <- simulateData(mean = betalm, SE = betalse)

#### Acquisition and inoculation estimates for High water stress
# Acquisition
alphahm <- 0.591; alphahse <- 0.058 # acquisition; probability; mean & SE -- MC normal distr
# Generate simulated data
alphahVec <- simulateData(mean = alphahm, SE = alphahse)

# Inoculation
betahm <- 0.450; betahse <- 0.047  # inoculation; probability; mean & SE -- MC normal distr
# Generate simulated data
betahVec <- simulateData(mean = betahm, SE = betahse)

#### Acquisition and inoculation estimates for extreme/severe water stress
# Acquisition
alphaexm <- 0.771; alphaexse <- 0.071 # acquisition; probability; mean & SE -- MC normal distr
# Generate simulated data
alphaexVec <- simulateData(mean = alphaexm, SE = alphaexse)

# Inoculation
betaexm <- 0.400; betaexse <- 0.219  # inoculation; probability; mean & SE -- MC normal distr
# Generate simulated data
betaexVec <- simulateData(mean = betaexm, SE = betaexse)


#### Incubation and latent period estimates are generated differently
#### Generate data sets with 1 column for each stress scenario
# Symptom onset rate (incubation rate, delta)
# Analyzing symptom data from choice experiment
choiceData <- read.csv("may1file.csv", header = TRUE)
summary(choiceData)
# remove outliers
choiceData <- choiceData[!choiceData$Populations == max(choiceData$Populations),]
choiceData <- choiceData[!choiceData$Populations == 0,]
choiceData$source.mpa <- choiceData$source.mpa/10 # Convert source plant water potential to -mPa
choiceData$treatment <- factor(choiceData$treatment, levels(choiceData$treatment)[c(3,2,1)])
# regression of symptoms and water treatment as an ordinal variable
sympModel <- glm(symptoms ~ treat*variety, data = choiceData, family = "poisson")
#plot(sympModel)
summary(sympModel)
predSymp <- predict(sympModel, type = "response", se.fit = TRUE)
choiceData$predSymp <- predSymp[[1]]
choiceData$predSE <- predSymp[[2]]
sympMeans <- aggregate(choiceData$predSymp, by = list(choiceData$treatment), mean, na.rm = TRUE)
sympErrors <- aggregate(choiceData$predSE, by = list(choiceData$treatment), mean, na.rm = TRUE)
sympMeans$SE <- sympErrors[,2]
names(sympMeans) <- c("water.level", "mean", "SE")
# Generate simulated data
deltaVec <- as.data.frame(apply(sympMeans[,2:3], 1, function(x) simulateData(mean = x[1], SE = x[2])))
names(deltaVec) <- sympMeans$water.level


# Exploring p and gamma parameter ranges
par(mfrow = c(3,1), mar = c(2,3,1,1))

for(i in 1:length(names(deltaVec))){
  waterlevel.i <- names(deltaVec)[i]
  hist(deltaVec[,waterlevel.i], xlim = c(0,5), main = waterlevel.i)
}

par(mfrow = c(1,1))


##################################################################################
# Combine parameter estimates into separate matrices and run model for each scenario
# Low stress parameters
paramLS <- cbind(alphalVec, betalVec, deltaVec$well, gamma, pVec, mu, nu)
# High stress parameters
paramHS <- cbind(alphahVec, betahVec, deltaVec$mod, gamma, pVec, mu, nu)
# Extreme stress parameters
paramEXS <- cbind(alphaexVec, betaexVec, deltaVec$severe, gamma, pVec, mu, nu)


#################################################################################
#### Loading cluster run data
WSDataList <- readRDS("water_stress_all_data_2016-01-28.rds")
Pepid <- WSDataList$Pepid
WSdata <- WSDataList$WSdata
WSTEdata <- WSDataList$WSTEdata

##################################################################################
#### Run simulations for equilibrial values for each scenario
LSRun <- as.data.frame(rbindlist(apply(paramLS, 1, WSsimFunc1)))
LSRun$stress.level <- "low"

HSRun <- as.data.frame(rbindlist(apply(paramHS, 1, WSsimFunc1)))
HSRun$stress.level <- "high"

EXSRun <- as.data.frame(rbindlist(apply(paramEXS, 1, WSsimFunc1)))
EXSRun$stress.level <- "severe"


#### Summaries of three runs
WSdata <- rbind(LSRun[,-1], HSRun[,-1], EXSRun[,-1])
WSdata$stress.level <- factor(WSdata$stress.level)
#saveRDS(WSdata, file = "water_stress_model_outputs_2015-11-11.rds")

#WSdata <- readRDS("water_stress_model_outputs_2016-01-28.rds")
WSdata$I <- round(WSdata$I, digits = 0)
#### Calculate the probability of an epidemic (i.e., proportion of total runs with I* > 0)
Pepid <- aggregate(WSdata$I, by = list(WSdata$stress.level), 
                   function(x) length(which(x > 0))/nrow(paramLS))
names(Pepid) <- c("stress.level", "Prob")
Pepid$se <- sqrt((Pepid$Prob*(1-Pepid$Prob))/nrow(paramLS))
Pepid$stress.level <- factor(Pepid$stress.level, levels(Pepid$stress.level)[c(2,1,3)])
Pepid <- Pepid[order(levels(Pepid$stress.level)),]

tiff(filename = "water stress prob epidemic 2016-01-28.tif") 
  with(Pepid,
       xyplot(Prob ~ stress.level, groups = stress.level,
              ly = Prob - se, uy = Prob + se,
              scales = list(col = 1, alternating = 1, tck = c(1, 0), 
                            x = list(cex = 1.3, labels = rep("",3)),
                            y = list(limits = c(0.5, 1), cex = 1.3,
                                     at = seq(0.5, 1, by = 0.1))
              ),
              xlab = "", 
              ylab = list("Probability of an epidemic", cex = 1.5),
              layout = c(1,1), pch = c(16), cex = 1.5,
              type = 'p', aspect = 1, col = "black",
              prepanel = prepanel.ci,                      
              panel = function(x, y, ...) {                
                panel.abline(v = unique(as.numeric(x)),  
                             col = "white")              
                panel.superpose(x, y, ...)               
              },                                          
              panel.groups = panel.ci))
dev.off()


#### Calculate mean and CI for I and V
WSEpidData <- WSdata[WSdata$I > 0,-3]
#hist(WSEpidData$I, breaks = seq(0,200,10))
#meanC <- WSEpidData %>% group_by(stress.level) %>% summarise(.,mean = mean(C), cil = quantile(C,0.0275), ciu = quantile(C,0.975))
meanI <- WSEpidData %>% group_by(stress.level) %>% summarise(.,mean = mean(I), cil = quantile(I,0.0275), ciu = quantile(I,0.975))
meanV <- WSEpidData %>% group_by(stress.level) %>% summarise(.,mean = (mean(V)/M)*100, cil = (quantile(V,0.0275)/M)*100, ciu = (quantile(V,0.975)/M)*100)
meanSpread <- rbind(meanV,meanI)
meanSpread$State <- c(rep("V",3),rep("I",3))
meanSpread$stress.level <- factor(meanSpread$stress.level, levels(meanSpread$stress.level)[c(2,1,3)])
meanSpread$dummyx <- ifelse(meanSpread$State == "V", as.numeric(meanSpread$stress.level) - 0.1, as.numeric(meanSpread$stress.level) + 0.1)

tiff(filename = "water stress mean infected host density 2016-01-29.tif") 
  with(meanSpread,
       xyplot(mean ~ dummyx, groups = State,
              ly = cil, uy = ciu,
              scales = list(col = 1, alternating = 1, tck = c(1, 0), 
                            x = list(limits = c(0.5,3.5), cex = 1.3, 
                                     at = seq(1,3,by=1), labels = levels(meanSpread$stress.level)),
                            y = list(limits = c(0, 100), cex = 1.3,
                                     at = seq(0, 100, by = 10))
              ),
              xlab = list("Water stress level", cex = 1.5), 
              ylab = list("Percent infectious or diseased", cex = 1.5),
              layout = c(1,1), pch = c(16), cex = 1.5,
              type = 'p', aspect = 1, col = c("darkgrey", "black"),
              prepanel = prepanel.ci,                      
              panel = function(x, y, ...) {                
                panel.abline(v = unique(as.numeric(x)),  
                             col = "white")              
                panel.superpose(x, y, ...)               
              },                                          
              panel.groups = panel.ci))
dev.off()



#########################################################################################
#### Calculating time to equilibrium
# Run simulated data for each scenario
LSTime <- as.data.frame(t(apply(paramLS, 1, WSTimeEqFunc)))
LSTime$stress.level <- "low"

HSTime <- as.data.frame(t(apply(paramHS, 1, WSTimeEqFunc)))
HSTime$stress.level <- "high"

EXSTime <- as.data.frame(t(apply(paramEXS, 1, WSTimeEqFunc)))
EXSTime$stress.level <- "severe"

#### Summaries of three runs
WSTEdata <- rbind(LSTime, HSTime, EXSTime)
names(WSTEdata) <- c("Istar", "TEQI", "Vstar", "TEQV", "stress.level")
WSTEdata$stress.level <- factor(WSdata$stress.level)
#saveRDS(WSTEdata, file = "water_stress_model_time_to_equilibrium_2015-11-11.rds")

#WSTEdata <- readRDS("water_stress_model_time_to_equilibrium_2016-01-28.rds")
# Round to whole numbers/"individuals"
WSTEdata$Istar <- round(WSTEdata$Istar, digits = 0)
WSTEdata$Vstar <- round(WSTEdata$Vstar, digits = 0)
# Remove simulations where no epidemic occurred (I* = 0)
WSTEdata <- WSTEdata[WSTEdata$Istar > 0,]

## Summary statistics
meanTEQ <- summarise(group_by(WSTEdata, stress.level), median = median(TEQI), cil = quantile(TEQI, 0.0275), ciu = quantile(TEQI, 0.975))
meanTEQ$stress.level <- factor(meanTEQ$stress.level, levels(meanTEQ$stress.level)[c(2,1,3)])
meanTEQ <- meanTEQ[order(levels(meanTEQ$stress.level)),]

tiff(filename = "water stress time eq plot 2016-01-28-2.tif") 
  with(meanTEQ,
       xyplot(mean~stress.level, groups = stress.level,
              ly = cil, uy = ciu,
              scales = list(col = 1, alternating = 1, tck = c(1, 0), 
                            x = list(cex = 1.3, labels = levels(meanTEQ$stress.level)),
                            y = list(limits = c(0, 250), cex = 1.3,
                                     at = seq(0, 250, by = 50))
              ),
              xlab = list("Water stress level", cex = 1.5), 
              ylab = list("Time to equilibrium", cex = 1.5),
              layout = c(1,1), pch = c(16), cex = 1.5,
              type = 'p', aspect = 1, col = "black",
              prepanel = prepanel.ci,                      
              panel = function(x, y, ...) {                
                panel.abline(v = unique(as.numeric(x)),  
                             col = "white")              
                panel.superpose(x, y, ...)               
              },                                          
              panel.groups = panel.ci))
dev.off()


###############################################################################
#### Estimating R0 for each water stress level
paramList <- list(paramLS, paramHS, paramEXS)
names(paramList) <- c("low", "high", "severe")
R0list <- list()
for(i in 1:length(paramList)){
  stress.level.i <- names(paramList)[i]
  params.i <- paramList[[i]]
  R0.i <- apply(params.i, 1, R0func)
  R0list[[i]] <- data.frame(R0 = R0.i, stress.level = stress.level.i)
}
R0data <- as.data.frame(rbindlist(R0list))

meanR0 <- summarise(group_by(R0data, stress.level), median = median(R0), cil = quantile(R0, 0.0275), ciu = quantile(R0, 0.975))
meanR0$stress.level <- factor(meanR0$stress.level, levels(meanR0$stress.level)[c(2,1,3)])
meanR0 <- meanR0[order(levels(meanR0$stress.level)),]

