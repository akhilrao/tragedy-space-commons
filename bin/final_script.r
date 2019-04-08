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
#font_import(prompt=FALSE)

#############################################################################
# 1a. Run calibration scripts, enable JIT compilation, adjust affinity mask, load functions and algorithms
#############################################################################

find_best_nls_parms <- 0 # 1: grid search to find the best starting values for NLS
physics_bootstrap <- 0 # 1: run the physics sensitivity analysis
source("calibrate_physical_model.r")

setwd("../bin/")
source("calibrate_econ_model.r")

rm(list=ls())
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
# 1b. Read in command args
#############################################################################

args <- commandArgs(trailingOnly=TRUE)

#############################################################################
# 1c. Set computation hyperparameters
#############################################################################

upper <- 1e6 # upper limit for some rootfinders - only requirement is that it should never bind
ncores <- 30 # number of cores to use for parallel computations
oa_gridsize <- 35
S_gridsize_opt <- 35
D_gridsize_opt <- 35
S_grid_upper_oa <- 8000
S_grid_upper_opt <- 8000
D_grid_upper_oa <- 250000
D_grid_upper_opt <- 25000

D_fraction_to_remove <- 0.5 # fraction of debris removed every period once removal is online. have no removal, set D_fraction_to_remove to 0.
D_removal_start_year <- 2030#as.numeric(args[1]) # pick a year within the projection time frame.

bootstrap <- 0 # 1: run sensitivity analysis for tax path
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
source("calibrate_parameters.r")

#############################################################################
# 2. Compute sequences of open access and optimal policies
#############################################################################

R_frac <- 0 # leave off for baseline
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
for(rs_year in 2021:2034) {
	R_start_year <- rs_year
	if(removal_comparison==1){
		R_frac <- D_fraction_to_remove
		source("main_model_projection.r")
		source("removal_comparison_figures.r")
	}
}
#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

if(bootstrap == 1){
	source("main_model_bootstrap.r")
}
