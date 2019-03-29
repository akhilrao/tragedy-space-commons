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
library(extrafont)
font_import()

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

upper <- 1e15 # upper limit for some rootfinders - only requirement is that it should never bind
ncores <- 3 # number of cores to use for parallel computations
oa_gridsize <- 32
S_gridsize_opt <- 32
D_gridsize_opt <- 32
S_grid_upper_oa <- 8000
S_grid_upper_opt <- 8000
D_grid_upper_oa <- 250000
D_grid_upper_opt <- 25000

bootstrap <- 0 # 1: run sensitivity analysis for tax path
removal_comparison <- 0 # 1: compare baseline model to model with debris removal

total_time <- proc.time()[3]

#############################################################################
# 1d. Calibration
#############################################################################

# Setting the end_year different from the projection_end extends the Morgan Stanley revenue and total value projections an additional 5 years, to avoid any end-of-horizon effects for a forecast out to end_year (e.g. numerical distortions in steady-state value functions). The idea is to "project" out to projection_end using the mean annual growth rate of the Morgan Stanley projections, then truncate back to end_year to avoid any end-of-horizon effects.
start_year <- 2006 # beginning of simulation
end_year <- 2040 # final year for plots
projection_end <- 2050 # final year for calculation
source("calibrate_parameters.r")

R_frac <- 0.5 # fraction of debris removed every period once removal is online.
R_start_year <- 2028 # pick a year within the projection time frame. to turn it off, set R_frac to 0.
R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year)

#############################################################################
# 2. Compute sequences of open access and optimal policies
#############################################################################

source("main_model_estimation.r")

#############################################################################
# 3. Generate open access and optimal time paths
#############################################################################

#opt_start_year <- c(start_year,2010,2015,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035)
#opt_start_year <- c(start_year,2010,2015)
opt_start_year <- c(start_year,2010,2015,2020,2025,2030,2035)
source("main_model_projection.r")

#############################################################################
# 4. Draw plots, write output
#############################################################################

source("main_model_figures.r")

if(removal_comparison==1){
	source("removal_comparison_figures.r")
}
#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

if(bootstrap == 1){
	source("main_model_bootstrap.r")
}