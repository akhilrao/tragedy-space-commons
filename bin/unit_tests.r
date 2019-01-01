##### Unit tests for this project

#####
# Testing the optimal policy solver
#####
rm(list=ls())
library(pracma)
library(data.table)
library(rootSolve)
library(gridExtra)
library(ggplot2)
library(viridis)
library(doParallel)
library(progress)
library(plot3D)
library(reshape2)
library(fields)
source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

# build grid, generate guesses, initialize dynamic_vfi_solver output list
gridsize <- 48
gridlist <- build_grid(gridmin=0, gridmax=30000, gridsize, cheby=1)
### shift grid down to zero if chebysheving it moved it up - should be unnecessary with expanded chebysehv array (secant factor)
if(min(gridlist$base_piece)>0) {
	gridlist$base_piece <- gridlist$base_piece - min(gridlist$base_piece)
	gridlist$igrid <- gridlist$igrid - min(gridlist$igrid)
}
dvs_output <- list()

# define T, sequence of p, F, and asats
T <- 15
p <- rep(1,length=T)
#F <- c(rep(13,length=T/2),rep(10,length=T/2))
F <- seq(from=3,to=3,length.out=T)
#F <- rep(13,length=T)
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
launch_constraint <- rep(30000,length=T)
discount_rate <- 0.05#0.6660755
discount_fac <- 1/(1+discount_rate)

# Build initial guess of value function
S_T_1 <- gridlist$igrid$sats
D_T_1 <- gridlist$igrid$debs
S_T <- S_T_1*(1-L(S_T_1,D_T_1))
vguess <- matrix(S_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)

# Check that path solver works in all periods
dvs_output <- policy_function_path_solver(gridpanel,gridlist,asats,T,p,F,ncores=4)

check <- dvs_output[[5]]

# bind the list of solved policies into a long dataframe
policy_path <- rbindlist(dvs_output)

# compare tps time series to solved time series
grid_lookup <- data.frame(sats=policy_path$satellites,debs=policy_path$debris,F=policy_path$F)

tps_path <- tps_opt_path(0,0,p,F,policy_path,asats,launch_constraint,grid_lookup,ncores=4)
View(tps_path)
li_path <- linint_opt_path(0,0,p,F,policy_path$optimal_launch_pfn,asats,grid_lookup)
View(li_path)

tps_path_base <- ggplot(tps_path, aes(x=time)) + xlab("Time") + theme_minimal()
tps_launches <- tps_path_base + geom_line(aes(y=launches), size=1) + ylab("Launches")
tps_sats <- tps_path_base + geom_line(aes(y=satellites), size=1) + ylab("Satellites")
tps_debs <- tps_path_base + geom_line(aes(y=debris), size=1) + ylab("Debris")
tps_risk <- tps_path_base + geom_line(aes(y=collision_rate), size=1) + ylab("Collision risk")
tps_returns <- tps_path_base + geom_line(aes(y=returns), size=1) + ylab("Returns")
tps_costs <- tps_path_base + geom_line(aes(y=costs), size=1) + ylab("Costs")

grid.arrange(tps_launches,tps_sats,tps_debs,tps_risk,tps_returns,tps_costs,nrow=3,ncol=2)

png(file=paste0("timepath_",gridsize,"_pt_basegrid.png"))
grid.arrange(tps_launches,tps_sats,tps_debs,tps_risk,tps_returns,tps_costs,nrow=3,ncol=2)
dev.off()
