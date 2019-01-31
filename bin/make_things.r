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
#library(compiler)
source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.
#enableJIT(3) # turn on JIT compilation for all functions

#############################################################################
# Calibration
#############################################################################

risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
implied_econ_series <- read.csv("../data/implied_costs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

ncores=16
T <- nrow(econ_series)
upper <- 1e15
#upper_seq <- seq(from=1,by=1,length=T)

start_year <- 2006
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
asats <- observed_time_series$num_destr_asat[which(observed_time_series$year>=start_year)]

discount_rate <- 0.05#0.6660755
discount_fac <- 1/(1+discount_rate)

fe_eqm <- econ_coefs[1,2] + econ_coefs[2,2]*econ_series$r_s + econ_coefs[3,2]*econ_series$Ft_Ft
fe_eqm <- observed_time_series$risk[which(observed_time_series$year>2005)]
#p <- rep(1,length=length(fe_eqm))

comb_econ_series <- merge(econ_series,implied_econ_series,by=c("year"))

# take the raw observed cost and revenue values
# physecon <- merge(observed_time_series,econ_series,by=c("year","risk"))
# p <- physecon$tot_rev/physecon$payloads_in_orbit
# F <- physecon$tot_cost/physecon$payloads_in_orbit

# take the theory-adjusted cost and revenue values. The final period gets truncated because the theory-adjustment formula uses up one year, so pad the final observation with the last value in the series. For p, this value is observed. For F, this value is just a copy of the newly-penultimate value.
physecon <- merge(observed_time_series,comb_econ_series,by=c("year","risk"))
p <- physecon$pi_t/physecon$payloads_in_orbit
p <- c(p,econ_series$tot_rev[length(econ_series$tot_rev)]/observed_time_series$payloads_in_orbit[length(observed_time_series$payloads_in_orbit)])
F <- physecon$F_hat/physecon$payloads_in_orbit
F <- c(F,F[length(F)])

observed_launches <- observed_time_series$launch_successes[which(observed_time_series$year>=start_year)]+observed_time_series$launch_failures[which(observed_time_series$year>=start_year)]
launch_constraint <- cummax(observed_launches)

#make_things <- function(opt_gridsize,...) {
opt_gridsize <- 16

#############################################################################
# Open access time path
#############################################################################

# build grid, generate guesses, initialize solver output list
gridsize <- 32
gridlist <- build_grid(gridmin=0, gridmax=25000, gridsize, cheby=1)

S_T_1 <- gridlist$igrid$sats
D_T_1 <- gridlist$igrid$debs
S_T <- S_T_1*(1-L(S_T_1,D_T_1))
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
oa_dvs_output <- list()

# Run path solver
oa_dvs_output <- oa_pvfn_path_solver(oa_dvs_output,gridpanel,gridlist,asats,T,p,F,fe_eqm,ncores=ncores)

base_grid <- unique(gridlist$igrid[,1])
oa_vfn <- matrix(oa_dvs_output[[T]]$oa_fleet_vfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)	
oa_pfn <- matrix(oa_dvs_output[[T]]$oa_launch_pfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
dev.new(width=8,height=7,unit="in")
png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_oa_vfn_pfn",".png"))
par(mfrow=c(2,2))
image2D(z=oa_vfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("value function"))
image2D(z=oa_vfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("value function"))
image2D(z=oa_pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("policy function"))
image2D(z=oa_pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("policy function"))
dev.off()

# bind the list of solved policies into a long dataframe
pvfn_path <- rbindlist(oa_dvs_output)

# compare tps time series to solved time series
grid_lookup <- data.frame(sats=pvfn_path$satellites,debs=pvfn_path$debris,F=pvfn_path$F)

tps_path <- tps_path_gen(S0,D0,p,F,pvfn_path,asats,launch_constraint,grid_lookup,ncores=ncores,OPT=0)

oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(tps_path)),tps_path)

##### compare fit of generated series against actual
oa_merged <- merge(oa_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

oa_merged <- cbind(oa_merged,fe_eqm=fe_eqm)

oa_fitcomp_base <- ggplot(data=oa_merged,aes(x=year))
oa_launch_fit <- oa_fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
oa_sat_fit <- oa_fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
oa_deb_fit <- oa_fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
oa_risk_fit <- oa_fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

dev.new()
#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(oa_launch_fit,oa_sat_fit,oa_risk_fit,oa_deb_fit,ncol=2)
#dev.off()
View(oa_merged)
#############################################################################
# Optimal time path
#############################################################################

# build grid, generate guesses, initialize dynamic_vfi_solver output list
gridsize <- opt_gridsize
gridlist <- build_grid(gridmin=0, gridmax=25000, gridsize, cheby=1)

S_T_1 <- gridlist$igrid$sats
D_T_1 <- gridlist$igrid$debs
S_T <- S_T_1*(1-L(S_T_1,D_T_1))
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
opt_dvs_output <- list()

# Run path solver
opt_dvs_output <- opt_pvfn_path_solver(opt_dvs_output,gridpanel,gridsize,gridlist,asats,T,p,F,ncores=ncores)

base_grid <- unique(gridlist$igrid[,1])
print(base_grid)
opt_vfn <- matrix(opt_dvs_output[[T]]$optimal_fleet_vfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)	
opt_pfn <- matrix(opt_dvs_output[[T]]$optimal_launch_pfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)

dev.new(width=8,height=7,unit="in")
png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_vfn_pfn",".png"))
par(mfrow=c(2,2))
image2D(z=opt_vfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("value function"))
image2D(z=opt_vfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("value function"))
image2D(z=opt_pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("policy function"))
image2D(z=opt_pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("policy function"))
dev.off()

# bind the list of solved policies into a long dataframe
pvfn_path <- rbindlist(opt_dvs_output)

# compare tps time series to solved time series
grid_lookup <- data.frame(sats=pvfn_path$satellites,debs=pvfn_path$debris,F=pvfn_path$F)

tps_path <- tps_path_gen(S0,D0,p,F,pvfn_path,asats,launch_constraint,grid_lookup,ncores=ncores,OPT=1)

opt_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(tps_path)),tps_path)


##### compare fit of generated series against actual
opt_merged <- merge(opt_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

opt_merged <- cbind(opt_merged,fe_eqm=fe_eqm)

opt_fitcomp_base <- ggplot(data=opt_merged,aes(x=year))
opt_launch_fit <- opt_fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
opt_sat_fit <- opt_fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
opt_deb_fit <- opt_fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
opt_risk_fit <- opt_fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

dev.new()
#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(opt_launch_fit,opt_sat_fit,opt_risk_fit,opt_deb_fit,ncol=2)
#dev.off()
View(opt_merged)
#############################################################################
# Compare OA and OPT outcomes
#############################################################################

OA_OPT <- merge(oa_merged,opt_merged,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT <- merge(OA_OPT,econ_series,by=c("year"))

dev.new()
OA_OPT_base <- ggplot(data=OA_OPT,aes(x=year))
OA_OPT_launch <- OA_OPT_base + geom_line(aes(y=launches.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launches.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=launch_successes.oa),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (OPT:blue, OA:red) vs. observed (black)")
OA_OPT_sat <- OA_OPT_base + geom_line(aes(y=satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=payloads_in_orbit.oa),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
OA_OPT_deb <- OA_OPT_base + geom_line(aes(y=debris.sim.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.sim.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=debris.obs.oa),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
OA_OPT_risk <- OA_OPT_base + geom_line(aes(y=collision_rate.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=collision_rate.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=risk.oa),size=1) +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_historical_series.png"))
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
dev.off()

OA_OPT$riskPoA <- OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt # 1 represents open access reaching optimal efficiency, larger numbers show deviations (inefficiency)
OA_OPT$PoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa # 1 represents open access reaching optimal welfare, larger numbers show deviations (suboptimal oa welfare)
OA_OPT$flow_welfare_loss <- OA_OPT$fleet_flowv.oa - OA_OPT$fleet_flowv.opt
OA_OPT$npv_welfare_loss <- OA_OPT$fleet_vfn_path.oa - OA_OPT$fleet_vfn_path.opt
OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa - OA_OPT$collision_rate.opt)*F*1e+9

oaoptcomp_base <- ggplot(data=OA_OPT,aes(x=year))

risk_comps <- oaoptcomp_base + geom_line(aes(y=riskPoA),size=0.85) +
							geom_hline(yintercept=1,linetype="dashed",color="blue") +
							ylab("Price of Anarchy for launch decisions\n(1 is optimal, larger numbers = worse)") + xlab("year") + theme_minimal()

flow_welf_loss <- oaoptcomp_base + geom_line(aes(y=flow_welfare_loss),size=0.85) +
							geom_hline(yintercept=0,linetype="dashed",color="blue") +
							ylab("Flow welfare loss from open access\n(undiscounted $1b)") + xlab("year") + theme_minimal() +
							ggtitle("Historical cost of open access and optimal tax path")

npv_welf_loss <- oaoptcomp_base + geom_line(aes(y=npv_welfare_loss),size=0.85) +
							ylab("NPV welfare loss from open access ($1b)") + xlab("year") + theme_minimal() +
							ggtitle("")

opt_tax_path <- oaoptcomp_base + geom_line(aes(y=opt_tax_path),size=0.85) +
							ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_minimal()

dev.new()
grid.arrange(flow_welf_loss,npv_welf_loss,risk_comps,opt_tax_path,ncol=2)

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_historical_cost_tax.png"))
grid.arrange(flow_welf_loss,npv_welf_loss,risk_comps,opt_tax_path,ncol=2)
dev.off()

#}

# make_things(opt_gridsize=8)
# make_things(opt_gridsize=16)
# make_things(opt_gridsize=32)
# make_things(opt_gridsize=64)
# make_things(opt_gridsize=128)
# make_things(opt_gridsize=256)
