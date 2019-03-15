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

#############################################################################
# 1a. Run calibration scripts, enable JIT compilation, adjust affinity mask, load functions and algorithms
#############################################################################

find_best_nls_parms <- 0
physics_bootstrap <- 0
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

upper <- 1e15 # upper limit for some rootfinders - should never bind
ncores <- 3#as.numeric(args[1]) # number of cores to use for parallel computations
oa_gridsize <- 24
# what's the right size?
S_gridsize_opt <- 24#as.numeric(args[2]) 
D_gridsize_opt <- 24#as.numeric(args[2]) 
S_grid_upper_oa <- 8000
S_grid_upper_opt <- 8000
D_grid_upper_oa <- 250000
D_grid_upper_opt <- 25000

bootstrap <- 0

total_time <- proc.time()[3]

#############################################################################
# 1d. Calibration
#############################################################################

source("calibrate_parameters.r")

R_frac <- 0 # fraction of debris removed every period once removal is online.
R_start_year <- 2025 # pick a year within the projection time frame. to turn it off, set R_frac to 0.
R_start <- which(seq(from=start_year,by=1,length.out=T)==R_start_year)

#############################################################################
# 2. Compute sequences of open access and optimal policies
#############################################################################

source("main_model_estimation.r")

#############################################################################
# 3. Generate open access and optimal time paths
#############################################################################

opt_start_year <- c(start_year,2010,2015,2020,2025,2030,2035)
source("main_model_projection.r")

#############################################################################
# 4. Draw plots, write output
#############################################################################

source("main_model_figures.r")

#############################################################################
# 5. Load bootstrapped data and generate projection uncertainty plot
#############################################################################

if(bootstrap == 1){
	source("main_model_bootstrap.r")
}