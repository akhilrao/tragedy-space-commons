##### Script to generate main model results for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Compute sequences of open access policies
# 2. Compute sequences of optimal policies

#############################################################################
# 1a. Optimal policies and values
#############################################################################

opt_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_opt, Dgridmax=D_grid_upper_opt, Sgridlength=S_gridsize_opt, Dgridlength=D_gridsize_opt, cheby=1)

# generate value and policy guesses - use terminal period. keep this separate from the open access guesses to allow for different gridsizes.
S_T_1 <- opt_gridlist$igrid$sats
D_T_1 <- opt_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay
#V_T <- p[T]*S_T
V_T <- (p[T]*S_T_1 - F[T]*L(S_T_1,D_T_1)*S_T_1)/(1-discount_fac)
vguess <- matrix(V_T,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
lpguess <- matrix(0,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
gridpanel <- grid_to_panel(opt_gridlist,lpguess,vguess)

# initialize solver output list
opt_dvs_output <- list()

# run path solver
optimal_model_compute_time <- proc.time()[3]
message("Generating optimal models...")
sink("log.solve.txt", append=FALSE)
opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,S_gridsize_opt,D_gridsize_opt,opt_gridlist,asats,T,p,F,ncores=ncores))
sink()
optimal_model_compute_time <- round(proc.time()[3]-optimal_model_compute_time,3)
message("Time to generate optimal models: ",optimal_model_compute_time," seconds")

# bind the list of solved policies into a long dataframe
opt_pvfn_path <- rbindlist(opt_dvs_output)

# Diagnostics: view solved launch policy functions
# function(vfn,launch_pfn,Sbasegrid,Dbasegrid,labels)
kk <- 44
year_being_examined <- 2006+kk
vfn <- opt_dvs_output[[kk]]$opt_fleet_vfn
pfn <- opt_dvs_output[[kk]]$opt_launch_pfn
plot_pfn_vfn(vfn,pfn,opt_gridlist$S_base_piece,opt_gridlist$D_base_piece,labels=c(paste0("Value function in ",year_being_examined),paste0("Policy function in ",year_being_examined)))

# # Diagnostics: view spline interpolation
# current_sats <- opt_dvs_output[[kk]]$satellites
# current_debs <- opt_dvs_output[[kk]]$debris
# tps_x <- as.matrix(cbind(current_sats,current_debs))
# tps_y <- as.matrix(pfn)
# #tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y,lambda=0))
# tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y))
# spline_pfn_int <- as.vector(predict(tps_model,x=tps_x))
# plot_pfn_vfn(spline_pfn_int,pfn,opt_gridlist$S_base_piece,opt_gridlist$D_base_piece,labels=c(paste0("Smoothed policy function in ",year_being_examined),paste0("Raw policy function in ",year_being_examined)))

#############################################################################
# 1b. Open access policies and values
#############################################################################

# build grid
gridsize <- oa_gridsize
oa_gridlist <- build_grid(gridmin=0, Sgridmax=S_grid_upper_oa, Dgridmax=D_grid_upper_oa, Sgridlength=gridsize, Dgridlength=gridsize, cheby=1)

# generate value and policy guesses - use final period
S_T_1 <- oa_gridlist$igrid$sats
D_T_1 <- oa_gridlist$igrid$debs
S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay 
V_T <- p[T]*S_T
vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(oa_gridlist,lpguess,vguess)

# initialize solver output list
oa_dvs_output <- list()

# run path solver
openaccess_model_compute_time <- proc.time()[3]
message("Generating open-access models...")
sink("log.solve.txt", append=TRUE)
oa_dvs_output <- suppressWarnings(oa_pvfn_path_solver(oa_dvs_output,gridpanel,oa_gridlist,asats,T,p,F,fe_eqm,ncores=ncores))
sink()
openaccess_model_compute_time <- round(proc.time()[3]-openaccess_model_compute_time,3)
message("Time to generate open-access models: ",openaccess_model_compute_time," seconds")

# bind the list of solved policies into a long dataframe
oa_pvfn_path <- rbindlist(oa_dvs_output)
