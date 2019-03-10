##### Script to generate results for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Run calibration scripts, read in hyperparameters, load calibrated parameters
# 2. Compute sequences of open access and optimal policies
# 3. Generate time paths from policy sequences
# 4. Draw figures, write output

#############################################################################
# 0. Load packages
#############################################################################

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
library(compiler)

#############################################################################
# 1a. Run calibration scripts, enable JIT compilation, adjust affinity mask, load functions and algorithms
#############################################################################

find_best_nls_parms <- 0
physics_bootstrap <- 1
source("calibrate_physical_model.r")

setwd("../bin/")
source("calibrate_econ_model.r")

rm(list=ls())
enableJIT(3) # turn on JIT compilation for all functions
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")

#############################################################################
# 1b. Read in command args
#############################################################################

args <- commandArgs(trailingOnly=TRUE)

#############################################################################
# 1c. Set computation hyperparameters
#############################################################################

upper <- 1e15 # upper limit for some rootfinders - should never bind
ncores <- 3#as.numeric(args[1]) # number of cores to use for parallel computations
oa_gridsize <- 24
# what's the right size?
S_gridsize_opt <- 24#as.numeric(args[2]) 
D_gridsize_opt <- 24#as.numeric(args[2]) 
S_grid_upper_oa <- 8000
S_grid_upper_opt <- 8000
D_grid_upper_oa <- 500000
D_grid_upper_opt <- 50000

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  quiet(force(x)) 
} 

total_time <- proc.time()[3]

#############################################################################
# 1d. Calibration
#############################################################################

source("calibrate_parameters.r")

#############################################################################
# 2b. Optimal policies and values
#############################################################################

opt_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_opt, Dgridmax=D_grid_upper_opt, Sgridlength=S_gridsize_opt, Dgridlength=D_gridsize_opt, cheby=1) # gridmax=25000 seems to work well for the data

# generate value and policy guesses - use terminal period. keep this separate from the open access guesses to allow for different gridsizes.
S_T_1 <- opt_gridlist$igrid$sats
D_T_1 <- opt_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay #guess for BW parameterization
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
lpguess <- matrix(0,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
gridpanel <- grid_to_panel(opt_gridlist,lpguess,vguess)

# dev.new()
# par(mfrow=c(1,2))
# plot_pfn_vfn(vguess,lpguess,opt_gridlist$S_base_piece,opt_gridlist$D_base_piece,c("vfn","pfn"))

print(paste0("Satellite stock nodes: ",paste(c(round(opt_gridlist$S_base_piece,digits=5)), sep=" ")))
print(paste0("Debris stock nodes: ",paste(c(round(opt_gridlist$D_base_piece,digits=5)), sep=" ")))

# initialize solver output list
opt_dvs_output <- list()

# run path solver
print("Generating full-sample optimal models...")
sink("log.solve.txt", append=FALSE)
opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,S_gridsize_opt,D_gridsize_opt,opt_gridlist,asats,T,p,F,ncores=ncores))
sink()

# bind the list of solved policies into a long dataframe
opt_pvfn_path <- rbindlist(opt_dvs_output)

#############################################################################
# 2a. Open access policies and values
#############################################################################

# build grid
gridsize <- oa_gridsize
oa_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_oa, Dgridmax=D_grid_upper_oa, Sgridlength=gridsize, Dgridlength=gridsize, cheby=1)

# generate value and policy guesses - use final period
S_T_1 <- oa_gridlist$igrid$sats
D_T_1 <- oa_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay #guess for BW parameterization
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(oa_gridlist,lpguess,vguess)

print(paste0("Satellite stock nodes: ",paste(c(round(oa_gridlist$S_base_piece,digits=5)), sep=" ")))
print(paste0("Debris stock nodes: ",paste(c(round(oa_gridlist$D_base_piece,digits=5)), sep=" ")))

# initialize solver output list
oa_dvs_output <- list()

# run path solver
print("Calculating full-sample open access models...")
sink("log.solve.txt", append=TRUE)
oa_dvs_output <- suppressWarnings(oa_pvfn_path_solver(oa_dvs_output,gridpanel,oa_gridlist,asats,T,p,F,fe_eqm,ncores=ncores))
sink()

# bind the list of solved policies into a long dataframe
oa_pvfn_path <- rbindlist(oa_dvs_output)

#############################################################################
# 3. Generate open access and optimal time paths
#############################################################################

# specific_t <- T-2
# oa_specific_vfn <- matrix(oa_dvs_output[[specific_t]]$oa_fleet_vfn,nrow=oa_gridsize,ncol=oa_gridsize)
# oa_specific_pfn <- matrix(oa_dvs_output[[specific_t]]$oa_launch_pfn,nrow=oa_gridsize,ncol=oa_gridsize)
# opt_specific_vfn <- matrix(opt_dvs_output[[specific_t]]$opt_fleet_vfn,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
# opt_specific_pfn <- matrix(opt_dvs_output[[specific_t]]$opt_launch_pfn,nrow=S_gridsize_opt,ncol=D_gridsize_opt)

# dev.new()
# par(mfrow=c(2,2))
# plot_pfn_vfn(oa_specific_vfn,oa_specific_pfn,oa_gridlist$S_base_piece,oa_gridlist$D_base_piece,c("oa vfn","oa pfn"))
# plot_pfn_vfn(opt_specific_vfn,opt_specific_pfn,opt_gridlist$S_base_piece,opt_gridlist$D_base_piece,c("opt vfn","opt pfn"))

source("simulation_algorithms.r")
### open access paths
oa_grid_lookup <- data.frame(sats=oa_pvfn_path$satellites,debs=oa_pvfn_path$debris,F=oa_pvfn_path$F)
print("Generating full-sample open access path...")
oa_tps_path <- tps_path_gen(S0,D0,0,p,F,oa_pvfn_path,asats,launch_constraint,oa_grid_lookup,ncores=ncores,OPT=0,linear_policy_interp=0)
oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(oa_tps_path)),oa_tps_path)

### optimal paths
opt_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
print("Generating full-sample optimal path...")
opt_tps_path <- tps_path_gen(S0,D0,0,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
opt_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(opt_tps_path)),opt_tps_path)

#############################################################################
# 4. Draw plots, write output
#############################################################################

OA_OPT_full <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT_full <- merge(OA_OPT_full,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"),all=TRUE)

OA_OPT_full <- merge(OA_OPT_full,econ_series,by=c("year"),all=TRUE)

selected_years <- intersect(which(OA_OPT_full$year>=start_year),which(OA_OPT_full$year<=end_year))
OA_OPT <- OA_OPT_full[selected_years,]

# dev.new()
OA_OPT_base <- ggplot(data=OA_OPT[which(OA_OPT$year<=2015),],aes(x=year))
OA_OPT_launch <- OA_OPT_base + geom_line(aes(y=launches.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launches.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=launch_successes),size=1) +				
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated historical series \n(OPT:blue, OA:red) vs. observed (black)")
OA_OPT_sat <- OA_OPT_base + geom_line(aes(y=satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle(" \n ")
OA_OPT_deb <- OA_OPT_base + geom_line(aes(y=debris.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=debris),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
OA_OPT_risk <- OA_OPT_base + geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=risk.x/payloads_in_orbit),size=1) +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()
# grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
# dev.off()

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_historical_series.png"))
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
dev.off()


# dev.new()
OA_OPT_base <- ggplot(data=OA_OPT,aes(x=year))
OA_OPT_launch <- OA_OPT_base + geom_line(aes(y=launches.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launches.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=launch_successes),size=1) +					
							geom_vline(xintercept=2015,size=1,linetype="dashed") +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated historical and projected series\n(OPT:blue, OA:red) vs. observed (black)")
OA_OPT_sat <- OA_OPT_base + geom_line(aes(y=satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							geom_vline(xintercept=2015,size=1,linetype="dashed") +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle(" \n ")
OA_OPT_deb <- OA_OPT_base + geom_line(aes(y=debris.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=debris),size=1) +
							geom_vline(xintercept=2015,size=1,linetype="dashed") +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
OA_OPT_risk <- OA_OPT_base + geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color="red",size=0.8) +
							geom_line(aes(y=risk.x/payloads_in_orbit),size=1) +
							geom_vline(xintercept=2015,size=1,linetype="dashed") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()
# grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
# dev.off()

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_projected_series.png"))
grid.arrange(OA_OPT_launch,OA_OPT_sat,OA_OPT_risk,OA_OPT_deb,ncol=2)
dev.off()

# Price of Anarchy in terms of collision risk. 1 represents no loss to anarchy, larger numbers show larger losses from anarchy.
OA_OPT$riskPoA <- (OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt)*(OA_OPT$satellites.opt/OA_OPT$satellites.oa)
# Price of Anarchy in terms of flow welfare. 1 : no present gains or losses to anarchy, >1 : present losses to anarchy, <1 : present gains to anarchy.
OA_OPT$flowWelfPoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa 
# Price of Anarchy in terms of NPV of welfare. 1 : no permanent gains or losses to anarchy, >1 : permanent losses to anarchy, <1 : permanent gains to anarchy.
OA_OPT$NPVPoA <- OA_OPT$fleet_vfn_path.opt/OA_OPT$fleet_vfn_path.oa 

# Since we're using aggregate data we need to divide by the number of satellites to get things into per-satellite units. Dividing by open access # of satellites puts everything relative to the "initial condition of open access".
OA_OPT$flow_welfare_loss <- (OA_OPT$fleet_flowv.oa - OA_OPT$fleet_flowv.opt)*norm_const/OA_OPT$satellites.oa
OA_OPT$npv_welfare_loss <- (OA_OPT$fleet_vfn_path.oa - OA_OPT$fleet_vfn_path.opt)*norm_const/OA_OPT$satellites.oa

F_over_horizon <- F[1:nrow(OA_OPT)]
OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa/OA_OPT$satellites.oa - OA_OPT$collision_rate.opt/OA_OPT$satellites.opt)*F_over_horizon*norm_const*1e+9/OA_OPT$satellites.oa # 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

oaoptcomp_base <- ggplot(data=OA_OPT,aes(x=year))

risk_comps <- oaoptcomp_base + geom_line(aes(y=riskPoA),size=0.85) +
							geom_hline(yintercept=1,linetype="dashed",color="blue") +
							ylab("Price of Anarchy for launch decisions\n(1 is optimal, larger numbers = worse)") + xlab("year") + theme_minimal()

flow_welf_loss <- oaoptcomp_base + geom_line(aes(y=flow_welfare_loss),size=0.85) +
							geom_hline(yintercept=0,linetype="dashed",color="blue") +
							ylab("Flow welfare gap (open access-optimal) \n(undiscounted $1b)") + xlab("year") + theme_minimal() +
							ggtitle("Historical cost of open access and optimal tax path")

npv_welf_loss <- oaoptcomp_base + geom_line(aes(y=npv_welfare_loss),size=0.85) +
							ylab("NPV welfare loss from open access ($1b)") + xlab("year") + theme_minimal() +
							ggtitle("")

opt_tax_path <- oaoptcomp_base + geom_line(aes(y=opt_tax_path),size=0.85) +
							ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_minimal()

npv_poa_path <- oaoptcomp_base + geom_line(aes(y=NPVPoA),size=0.85) +
							geom_hline(yintercept=1,linetype="dashed",color="blue") +
							ylab("Price of Anarchy\n(1 is no gain, larger numbers = larger gain)") + xlab("year") + theme_minimal() +
							ggtitle("Improvement in global satellite fleet NPV\nfrom optimal management")
# dev.new()
# grid.arrange(flow_welf_loss,npv_welf_loss,risk_comps,opt_tax_path,ncol=2)
# dev.off()

png(width=700,height=700,filename=paste0("../images/",gridsize,"_pt_opt_simulated_historical_cost_tax.png"))
grid.arrange(flow_welf_loss,npv_welf_loss,risk_comps,opt_tax_path,ncol=2)
dev.off()

# dev.new()
# npv_poa_path
# dev.off()

png(width=400,height=400,filename=paste0("../images/",gridsize,"_pt_NPV_PoA_path.png"))
npv_poa_path
dev.off()

cat(paste0("\n Done. Total script wall time: ",round(proc.time()[3] - total_time,3)/60," minutes"))

write.csv(OA_OPT, file="../data/main_simulation.csv")

#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

source("final_script_bootstrap.r") # defaults to 36x36 grid, 30 bootstrap samples

main_sim <- read.csv("../data/main_simulation.csv")
bootstrap_sims <- read.csv("../data/bootstrapped_simulation.csv")

main_small <- data.frame(year=main_sim$year, 
						launches.oa=main_sim$launches.oa,
						satellites.oa=main_sim$satellites.oa,
						debris.oa=main_sim$debris.oa,
						collision_rate.oa=main_sim$collision_rate.oa,
						launches.opt=main_sim$launches.opt,
						satellites.opt=main_sim$satellites.opt,
						debris.opt=main_sim$debris.opt,
						collision_rate.opt=main_sim$collision_rate.opt,
						bootstrap_draw=0)

bootstrap_small <- data.frame(year=bootstrap_sims$year, 
						launches.oa=bootstrap_sims$launches.oa,
						satellites.oa=bootstrap_sims$satellites.oa,
						debris.oa=bootstrap_sims$debris.oa,
						collision_rate.oa=bootstrap_sims$collision_rate.oa,
						launches.opt=bootstrap_sims$launches.opt,
						satellites.opt=bootstrap_sims$satellites.opt,
						debris.opt=bootstrap_sims$debris.opt,
						collision_rate.opt=bootstrap_sims$collision_rate.opt,
						bootstrap_draw=bootstrap_sims$bootstrap_draw)						

m_small_long <- reshape(main_small, idvar=c("year","bootstrap_draw"), times=(colnames(main_small)[-c(1,ncol(main_small))]), direction="long", v.names="m.value", varying=(colnames(main_small)[-c(1,ncol(main_small))]))
bs_small_long <- reshape(bootstrap_small, idvar=c("year","bootstrap_draw"), times=(colnames(bootstrap_small)[-c(1,ncol(bootstrap_small))]), direction="long", v.names="bs.value", varying=(colnames(bootstrap_small)[-c(1,ncol(bootstrap_small))]))
rownames(m_small_long) <- NULL
rownames(bs_small_long) <- NULL

m_bs_small_long <- merge(m_small_long, bs_small_long, by=c("year","time","bootstrap_draw"), suffixes=c(".main",".bs"), all=TRUE)
rownames(m_bs_small_long) <- NULL
head(m_bs_small_long)

m_bs_small_long_bootstrap_oalaunch_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="launches.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Open access launches") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("Residual bootstrap simulations,",B," draws\n(black: main model projections, light blue: bootstrap model projections)\nOA launch projections"))
m_bs_small_long_bootstrap_oasats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Open access satellites") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOA satellite projections"))
m_bs_small_long_bootstrap_oadebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Open access debris") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOA debris projections"))
m_bs_small_long_bootstrap_oacoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Open access collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOA collision rate projections"))

dev.new()
grid.arrange(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot, ncol=2)



m_bs_small_long_bootstrap_optlaunch_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="launches.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Optimal launches") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("Residual bootstrap simulations,",B," draws\n(black: main model projections, light blue: bootstrap model projections)\nOptimal launch projections"))
m_bs_small_long_bootstrap_optsats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Optimal satellites") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOptimal satellite projections"))
m_bs_small_long_bootstrap_optdebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Optimal debris") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOptimal debris projections"))
m_bs_small_long_bootstrap_optcoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Optimal collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw))) +
						theme_minimal() + ggtitle(paste0("\n\nOptimal collision rate projections"))

dev.new()
grid.arrange(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot, ncol=2)
