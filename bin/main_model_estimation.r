##### Script to generate main model results for "Tragedy of the Space Commons" paper. This script estimates the policy sequences given the calibrated parameters.
###
# Script flow:
# 1. Compute sequences of open access policies
# 2. Compute sequences of optimal policies

#############################################################################
# 1a. Optimal policies and values
#############################################################################

opt_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_opt, Dgridmax=D_grid_upper_opt, Sgridlength=S_gridsize_opt, Dgridlength=D_gridsize_opt, cheby=1) # gridmax=25000 seems to work well for the data

# generate value and policy guesses - use terminal period. keep this separate from the open access guesses to allow for different gridsizes.
S_T_1 <- opt_gridlist$igrid$sats
D_T_1 <- opt_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay # in the final period, the launch rate is zero
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
lpguess <- matrix(0,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
gridpanel <- grid_to_panel(opt_gridlist,lpguess,vguess)

# initialize solver output list
opt_dvs_output <- list()

# run path solver
message("\nCalculating optimal models...")
sink("log.solve.txt", append=FALSE)
opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,S_gridsize_opt,D_gridsize_opt,opt_gridlist,asats,T,p,F,ncores=ncores))
sink()

# bind the list of solved policies into a long dataframe
opt_pvfn_path <- rbindlist(opt_dvs_output)

#############################################################################
# 1b. Open access policies and values
#############################################################################

# build grid
gridsize <- oa_gridsize
oa_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_oa, Dgridmax=D_grid_upper_oa, Sgridlength=gridsize, Dgridlength=gridsize, cheby=1)

# generate value and policy guesses - use final period
S_T_1 <- oa_gridlist$igrid$sats
D_T_1 <- oa_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay # guess for final period value
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(oa_gridlist,lpguess,vguess)

# initialize solver output list
oa_dvs_output <- list()

# run path solver
message("\nCalculating open access models...")
sink("log.solve.txt", append=TRUE)
oa_dvs_output <- suppressWarnings(oa_pvfn_path_solver(oa_dvs_output,gridpanel,oa_gridlist,asats,T,p,F,fe_eqm,ncores=ncores))
sink()

# bind the list of solved policies into a long dataframe
oa_pvfn_path <- rbindlist(oa_dvs_output)
