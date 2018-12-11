##### Unit tests for tihs project

#####
# Testing the optimal policy solver
#####
rm(list=ls())
library(rootSolve)
library(gridExtra)
library(ggplot2)
library(viridis)
library(doParallel)
library(progress)
source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")

# build grid and generate guesses
gridsize <- 24
gridlist <- build_grid(gridmin=0, gridmax=10000, gridsize, 1)
vguess <- matrix(gridlist$igrid$sats,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)

# define T, sequence of p, F, and asats
T <- 10
p <- rep(1,length=T)
F <- rep(10,length=T)
asats <- rep(0,length=T)

# define physical parameters and discount rate+factor
risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

start_year <- 2005
S0 <- observed_time_series$payloads_in_orbit[which(observed_time_series$year==start_year)]
D0 <- observed_time_series$debris[which(observed_time_series$year==start_year)]
avg_sat_decay <- mean(observed_time_series$payloads_decayed[which(observed_time_series$year>=start_year)]/observed_time_series$payloads_in_orbit[which(observed_time_series$year>=start_year)])

aS <- risk_cal[1,2]
aD <- risk_cal[2,2]
aSS <- risk_cal[3,2]
aSD <- risk_cal[4,2]
aDD <- risk_cal[5,2]

aDDbDD <- deblom_cal[7,2]
bSS <- deblom_cal[5,2]
bSD <- deblom_cal[6,2]
d <- deblom_cal[2,2]
Z_coef <- deblom_cal[1,2]
m <- deblom_cal[3,2]
asat_coef <- deblom_cal[4,2]
#asats <- observed_time_series$num_destr_asat[which(observed_time_series$year>=start_year)]

discount_rate <- 0.05#0.6660755
discount_fac <- 1/(1+discount_rate)

# Check that objective function makes sense
k=50
fleet_preval(X=0,S=gridpanel$S[k],D=gridpanel$D[k],value_fn=gridpanel$V,asats=asats,t=T,p=p,F=F,igrid=gridlist$igrid)

# Check that solver works in terminal period
#registerDoParallel(cores=2)
test <- dynamic_vfi_solver(gridpanel,igrid=gridlist$igrid,asats,T,T,p,F) 
#stopCluster()
