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

ncores <- 30 # number of cores to use for parallel computations
find_best_nls_parms <- 0 # 1: grid search to find the best starting values for NLS. takes some time, so is set to 0 and starts from prior solve results by default.
physics_bootstrap <- 1 # 1: run the physics sensitivity analysis.
source("plotting_functions.r")
# these scripts calibrate the parameters of the physical and economic models, and write the parameters out to separate files. they do not construct the data series' necessary for the value function iteration.
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.
source("calibrate_physical_model.r")
source("calibrate_econ_model.r")

rm(list=ls()) # clear workspace again, now that the models are calibrated
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
# 1b. Read in command args --- 071719: DEPRECATED, SCRIPT DOES NOT ACCEPT COMMAND ARGS
#############################################################################

args <- commandArgs(trailingOnly=TRUE)

#############################################################################
# 1c. Set computation hyperparameters
#############################################################################

ncores <- 30 # number of cores to use for parallel computations
upper <- 1e6 # upper limit for some rootfinders - only requirement is that it should never bind
oa_gridsize <- 35
S_gridsize_opt <- 35
D_gridsize_opt <- 35
S_grid_upper_oa <- 8000
S_grid_upper_opt <- 8000
D_grid_upper_oa <- 250000
D_grid_upper_opt <- 25000

D_fraction_to_remove <- 0.5 # fraction of debris removed every period once removal is online. have no removal, set D_fraction_to_remove to 0.
D_removal_start_year <- 2022 # pick a year within the projection time frame.

bootstrap <- 1 # 1: run sensitivity analysis for tax path. set to 1 by default.
n_path_sim_bootstrap_draws <- 125 # number of bootstrap draws to use for open access and optimal path sensitivity analysis. only matters when bootstrap <- 1.

removal_comparison <- 1 # 1: compare baseline model to model with debris removal. will (re)generate paths with R_frac <- 0.

total_time <- proc.time()[3]

#############################################################################
# 1d. Calibration
#############################################################################

# Setting the end_year different from the projection_end extends the Morgan Stanley revenue and total value projections an additional 5 years, to avoid any end-of-horizon effects for a forecast out to end_year (e.g. numerical distortions in steady-state value functions). The idea is to "project" out to projection_end using the mean annual growth rate of the Morgan Stanley projections, then truncate back to end_year to avoid any end-of-horizon effects.
start_year <- 2006 # beginning of simulation
end_year <- 2040 # final year for plots
projection_end <- 2050 # final year for calculation
opt_start_year <- c(start_year,2010,2015,2020,2025,2030,2035)
source("calibrate_parameters.r") # reads in all calibrated parameter values, estimates the launch constraint, constructs the necessary data series, and generates main text figure 1.

#############################################################################
# 2. Compute sequences of open access and optimal policies
#############################################################################

R_frac <- 0.5 # leave off for baseline
R_start_year <- D_removal_start_year
R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year)
source("main_model_estimation.r")

#############################################################################
# 3. Generate open access and optimal time paths
#############################################################################

source("main_model_projection.r")

#############################################################################
# 4. Draw plots, write output
#############################################################################

source("main_model_figures.r")

cat(paste0("\n Done. Total wall time for main model: ",round(proc.time()[3] - total_time,3)/60," minutes"))

total_time <- proc.time()[3]

if(removal_comparison==1){
	for(rs_year in 2021:2034) {
		R_start_year <- rs_year
		R_frac <- D_fraction_to_remove
		source("main_model_projection.r")
		# if the appropriate no-removal file exists, continue to generate figures. else, generate the appropriate no-removal file.
		if(file.exists(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))==FALSE) {
			R_frac <- 0
			source("main_model_projection.r")
			R_frac <- D_fraction_to_remove
		}
		#source("removal_comparison_figures.r") # 071519: FIX THIS SCRIPT, IT'S CAUSING FAILS OF THE TYPE 
		# Error in data.frame(list(npv_welfare_gain.rem = c(2472.32389737803, 1547.56435677332,  :   arguments imply differing number of rows: 5, 0

	}
}

cat(paste0("\n Done. Total wall time for removal models: ",round(proc.time()[3] - total_time,3)/60," minutes"))

#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

if(bootstrap == 1){
	R_frac <- 0 # 0 disables debris removal for the bootstrapped models
	source("main_model_bootstrap.r")
}
