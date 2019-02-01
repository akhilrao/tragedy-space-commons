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
# Set computation hyperparameters
#############################################################################

upper <- 1e15 # upper limit for some rootfinders - basically should never bind
ncores <- 16 # number of cores to use for parallel computations
oa_gridsize <- 32
opt_gridsize <- 32

#############################################################################
# Calibration
#############################################################################

source("calibrate_parameters.r")

#############################################################################
# Open access policies and values
#############################################################################

# build grid

# make_gridpanel <- function(gridmin,gridmax,gridsize,cheby) {
# 	gridsize <- oa_gridsize
# 	gridlist <- build_grid(gridmin=gridmin, gridmax=gridmax, gridsize, cheby=1)
# 	# generate value and policy guesses - use final period
# 	S_T_1 <- gridlist$igrid$sats
# 	D_T_1 <- gridlist$igrid$debs
# 	S_T <- S_T_1*(1-L(S_T_1,D_T_1))
# 	V_T <- p[T]*S_T
# 	vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
# 	lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
# 	gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
# 	return(list(gridlist=gridlist,gridpanel=gridpanel))
# }

# gridlist_1 <- make_gridpanel(0,25000,gridsize,1)
# gridlist_2 <- make_gridpanel(0,5000,gridsize,1)
# gridpanel <- rbind(gridlist_1$gridpanel,gridlist_2$gridpanel)
# gridlist_1 <- gridlist_1$gridlist
# gridlist_2 <- gridlist_2$gridlist

# gridlist <- list(base_piece=c(gridlist_1$base_piece,gridlist_2$base_piece), igrid=cbind(gridlist_1$igrid,gridlist_2$igrid))

gridsize <- oa_gridsize
gridlist <- build_grid(gridmin=0, gridmax=20000, gridsize, cheby=1)

# generate value and policy guesses - use final period
S_T_1 <- gridlist$igrid$sats
D_T_1 <- gridlist$igrid$debs
S_T <- S_T_1*(1-L(S_T_1,D_T_1))
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)

print(round(gridlist$base_piece,digits=5))

# initialize solver output list
oa_dvs_output <- list()

# run path solver
oa_dvs_output <- suppressWarnings(oa_pvfn_path_solver(oa_dvs_output,gridpanel,gridlist,asats,T,p,F,fe_eqm,ncores=ncores))

# bind the list of solved policies into a long dataframe
oa_pvfn_path <- rbindlist(oa_dvs_output)

#############################################################################
# Optimal policies and values
#############################################################################

# build grid
gridsize <- opt_gridsize # not clear what the right choice is here - as high as possible? larger values require a desktop with multiple cores and lots of RAM.
# another issue: larger values make the (0,0) launch decision much larger than the OA decision. is this a numerical artefact in the solve? or a result of parameter values which limit Kessler possibilities? HOW TO TEST THIS? i can maybe rule out Kessler with fake parms?
gridlist <- build_grid(gridmin=0, gridmax=20000, gridsize, cheby=1) # gridmax=25000 seems to work well for the data

# generate value and policy guesses - use terminal period. keep this separate from the open access guesses to allow for different gridsizes.
S_T_1 <- gridlist$igrid$sats
D_T_1 <- gridlist$igrid$debs
S_T <- S_T_1*(1-L(S_T_1,D_T_1))
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)

print(round(gridlist$base_piece,digits=5))

# initialize solver output list
opt_dvs_output <- list()

# run path solver
opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,gridsize,gridlist,asats,T,p,F,ncores=ncores))

# bind the list of solved policies into a long dataframe
opt_pvfn_path <- rbindlist(opt_dvs_output)

#############################################################################
# Generate open access and optimal paths
#############################################################################

### open access paths
oa_grid_lookup <- data.frame(sats=oa_pvfn_path$satellites,debs=oa_pvfn_path$debris,F=oa_pvfn_path$F)
oa_tps_path <- tps_path_gen(S0,D0,0,p,F,oa_pvfn_path,asats,launch_constraint,oa_grid_lookup,ncores=ncores,OPT=0,linear_policy_interp=0)
oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(oa_tps_path)),oa_tps_path)

### optimal paths
opt_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
opt_tps_path <- tps_path_gen(S0,D0,0,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
opt_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(opt_tps_path)),opt_tps_path)

#############################################################################
# Compare OA and OPT paths
#############################################################################

OA_OPT <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT <- merge(OA_OPT,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))
OA_OPT <- merge(OA_OPT,econ_series,by=c("year"))

dev.new()
OA_OPT_base <- ggplot(data=OA_OPT,aes(x=year))
OA_OPT_launch <- OA_OPT_base + geom_line(aes(y=launches.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launches.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (OPT:blue, OA:red) vs. observed (black)")
OA_OPT_sat <- OA_OPT_base + geom_line(aes(y=satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
OA_OPT_deb <- OA_OPT_base + geom_line(aes(y=debris.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=debris),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
OA_OPT_risk <- OA_OPT_base + geom_line(aes(y=collision_rate.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=collision_rate.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=risk.y),size=1) +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_historical_series.png"))
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
dev.off()

OA_OPT$riskPoA <- OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt # 1 represents open access reaching optimal efficiency, larger numbers show deviations (inefficiency)
OA_OPT$PoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa # 1 represents open access reaching optimal welfare, larger numbers show deviations (suboptimal oa welfare)
OA_OPT$flow_welfare_loss <- (OA_OPT$fleet_flowv.oa - OA_OPT$fleet_flowv.opt)*norm_const/OA_OPT$satellites.oa
OA_OPT$npv_welfare_loss <- (OA_OPT$fleet_vfn_path.oa - OA_OPT$fleet_vfn_path.opt)*norm_const/OA_OPT$satellites.oa
OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa - OA_OPT$collision_rate.opt)*F*1e+9*norm_const/OA_OPT$satellites.oa

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


