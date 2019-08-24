##### Script to generate results for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Run calibration scripts, read in hyperparameters, load calibrated parameters
# 2. Compute sequences of open access and optimal policies
# 3. Generate time paths from policy sequences
# 4. Draw figures and write output
# 5. Run bootstrap physics sensitivity analysis if desired

#############################################################################
# 0. Load packages
#############################################################################

rm(list=ls())

library(pracma)
library(data.table)
library(rootSolve)
library(grid)
library(gridExtra)
library(ggplot2)
library(viridis)
library(doParallel)
library(progress)
library(plot3D)
library(reshape2)
library(fields)
library(compiler)
library(stargazer)
library(cowplot)
library(tidyr)
library(plyr)
library(extrafont)
library(glmnet)
library(reshape2)
library(BB)
library(ggpubr)

#font_import(prompt=FALSE)

#############################################################################
# 1a. Run calibration scripts, enable JIT compilation, adjust affinity mask, load functions and algorithms
#############################################################################

ncores <- 3 # number of cores to use for parallel computations
find_best_nls_parms <- 0 # 1: grid search to find the best starting values for NLS. takes some time; default is set to 0 and starts from prior solve results.
physics_bootstrap <- 0 # 1: run the physical calibration sensitivity analysis again. only necessary if parameter sets are to be regenerated from scratch. takes some time; default is set to 0 and starts from prior solve results.
n_physical_bootstrap_draws <- 1000 # number of draws for physical calibration sensitivity analysis. default is 1000.

source("plotting_functions.r")
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

# these scripts estimate the parameters of the physical and economic models, and write the parameters out to csv files for calibration later on. they do not construct the data series' necessary for the value function iteration.
source("calibrate_physical_model.r")
source("calibrate_econ_model.r")

rm(list=ls()) # clear workspace again, now that the model parameters are estimated
enableJIT(3) # turn on JIT compilation for all functions
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  quiet(force(x)) 
} 

#############################################################################
# 1b. Set computation hyperparameters
#############################################################################

ncores <- 3 # number of cores to use for parallel computations
upper <- 1e6 # upper limit for some rootfinders - only requirement is that it should never bind
oa_gridsize <- 28#35
S_gridsize_opt <- 28#35
D_gridsize_opt <- 28#35
##### 072219: try shrinking these to remove the "jump at 2020" launch rate artifact
S_grid_upper_oa <- 8000 
S_grid_upper_opt <- 3000 #8000
D_grid_upper_oa <- 250000
D_grid_upper_opt <- 10000 #25000

bootstrap <- 0 # 1: run sensitivity analysis for model outputs. set to 1 by default.
n_path_sim_bootstrap_draws <- 50 # number of bootstrap draws to use for open access and optimal path sensitivity analysis. only matters when bootstrap <- 1.

removal_comparison <- 1 # 1: compare baseline model to model with debris removal. will generate paths with R_frac <- 0 if necessary.

total_time <- proc.time()[3]

#############################################################################
# 1c. Calibration
#############################################################################

# Setting the end_year different from the projection_end extends the Morgan Stanley revenue and total value projections an additional 5 years, to avoid any end-of-horizon effects for a forecast out to end_year (e.g. numerical distortions in steady-state value functions). The idea is to "project" out to projection_end using the mean annual growth rate of the Morgan Stanley projections, then truncate back to end_year to avoid any end-of-horizon numerical artifacts.
start_year <- 2006 # beginning of simulation
end_year <- 2040 # final year for plots
projection_end <- 2050 # final year for calculation. should be weakly greater than end_year.
opt_start_year <- c(start_year,2010,2015,2020,2025,2030,2035)
opt_start_year_bs <- c(2006,2020,2035) # optimal management start years for bootstrap draws. WARNING: each entry here will add a lot (n_path_sim_bootstrap_draws*(time to compute a single optimal model)) to runtime! expand list with caution! (or with abundant cheap compute.) default is 2020 and 2035, to generate histogram of npv_welfare_gains comparable to the headline numbers.
source("calibrate_parameters.r", print.eval=TRUE) # reads in all calibrated parameter values, estimates the launch constraint, constructs the necessary data series, and generates main text figure 1.

#############################################################################
# 2. Compute sequences of open access and optimal policies
#############################################################################

D_fraction_to_remove <- 0.5 # fraction of debris removed every period once removal is online. default is 0.5, for use inside removal_comparison loop. to have no removal, set to 0. this variable is the "master copy" which stays constant inside the removal_comparison loop. if removal_comparison==0, this is irrelevant.
R_frac <- 0 # fraction of debris removed every period once removal is online. default is 0, so that main_model_projection generates no-removal projection bootstraps. this variable gets updated with removal_comparison inner loops, and reset to the value of D_fraction_to_remove. NOTE: set this to >0 if you want main model projections/bootstraps with removal.

D_removal_start_year <- 2027 # pick a year within the projection time frame.
R_start_year <- D_removal_start_year # same deal as R_frac: this is the copy that gets updated in the removal_comparison inner loops.
R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year) # this gets the correct integer label for the chosen R_start_year, which is used in the projection algorithm 
source("main_model_estimation.r")

#############################################################################
# 3. Generate open access and optimal time paths
#############################################################################

source("main_model_projection.r")

#############################################################################
# 4. Draw plots, write output
#############################################################################

source("main_model_figures.r", print.eval=TRUE)

message(paste0("\n Done. Total wall time for main model: ",round((proc.time()[3] - total_time)/60,3)," minutes"))

total_time <- proc.time()[3]

if(removal_comparison==1){
	# this loop generates removal outcomes for all years between 2021:2034
	for(rs_year in 2021:2034) {
		R_start_year <- rs_year
		R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year)
		R_frac <- D_fraction_to_remove
		source("main_model_projection.r")
		# if the appropriate no-removal file exists, continue to generate figures. else, generate the appropriate no-removal file.
		if(file.exists(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))==FALSE) {
			R_frac <- 0
			source("main_model_projection.r")
			R_frac <- D_fraction_to_remove
		}
		source("removal_summary_stat_files.r") # script to generate the files necessary to compute debris removal summary statistics
	}
	# the next three lines generate the figures for the main removal model
	R_start_year <- D_removal_start_year
	R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year) 
	
	# this block calculates summary statistics for changes in fleet NPV caused by debris removal starting in different years
	coi_list <- list()
	for(coi_rs_year in 2021:2034) {
		coi_list[[(coi_rs_year-2020)]] <- read.csv(paste0("../data/2006_7_starts_remfrac_0.5_remstart_",coi_rs_year,"_coi_base_dfrm.csv"))
	}
	coi_total_dfrm <- rbindlist(coi_list)
	coi_total_summary <- data.frame(
		best_change=as.numeric(summary(coi_total_dfrm$coi_effect_of_removal_pc)[1]), 
		avg_change=as.numeric(summary(coi_total_dfrm$coi_effect_of_removal_pc)[4]), 
		worst_change=as.numeric(summary(coi_total_dfrm$coi_effect_of_removal_pc)[6]))

	write.csv(coi_total_summary, file="../data/pc_effect_of_removal_summary.csv")

	source("removal_comparison_figures.r", print.eval=TRUE)
}

message(paste0("\n Done. Total wall time for removal models: ",round((proc.time()[3] - total_time)/60,3)," minutes"))

#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

if(bootstrap == 1){
	R_frac <- 0 # 0 disables debris removal for the bootstrapped models
	source("main_model_bootstrap.r", print.eval=TRUE)
}
