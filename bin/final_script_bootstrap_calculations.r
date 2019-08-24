##### Script to generate results for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Run calibration scripts, read in hyperparameters, load calibrated parameters
# 2. Compute sequences of open access and optimal policies
# 3. Generate time paths from policy sequences
# 4. Draw figures, write output

#############################################################################
# 0. Begin loop
#############################################################################

# BEGIN BOOTSTRAP LOOP

for(b in 1:n_path_sim_bootstrap_draws) {

total_time <- proc.time()[3]

message(paste0("\n\nBeginning bootstrap draw ",b))
	#############################################################################
	# 1. Calibration
	#############################################################################

	source("calibrate_parameters_bootstrap.r")

	#############################################################################
	# 2a. Optimal policies and values
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
	opt_start_time <- proc.time()[3]
message(paste0("\nCalculating optimal models for draw ",b,"...\n"))
sink("log.solve.txt", append=TRUE)
	opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,S_gridsize_opt,D_gridsize_opt,opt_gridlist,asats,T,p,F,ncores=ncores))
sink()
	opt_end_time <- proc.time()[3]

	# bind the list of solved policies into a long dataframe
	opt_pvfn_path <- rbindlist(opt_dvs_output)

	#############################################################################
	# 2b. Open access policies and values
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
message(paste0("\nCalculating open access models for draw ",b,"...\n"))
sink("log.solve.txt", append=TRUE)
	oa_dvs_output <- suppressWarnings(oa_pvfn_path_solver(oa_dvs_output,gridpanel,oa_gridlist,asats,T,p,F,fe_eqm,ncores=ncores))
sink()

	# bind the list of solved policies into a long dataframe
	oa_pvfn_path <- rbindlist(oa_dvs_output)

	#############################################################################
	# 3. Generate open access and optimal time paths
	#############################################################################

	### open access paths
	oa_grid_lookup <- data.frame(sats=oa_pvfn_path$satellites,debs=oa_pvfn_path$debris,F=oa_pvfn_path$F)
message(paste0("\nGenerating open access path for draw ",b,"...\n"))
sink("log.solve.txt", append=TRUE)
	oa_tps_path <- tps_path_gen(S0,D0,0,R_start,R_start_year,R_frac,p,F,oa_pvfn_path,asats,launch_constraint,oa_grid_lookup,ncores=ncores,OPT=0,linear_policy_interp=0)
sink()
	oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(oa_tps_path)),oa_tps_path)

	### optimal paths -- 081319: FIX THIS, THE 2006 RUNS AREN'T BEING GENERATED
	opt_path_list <- list()
for(o in 1:length(opt_start_year_bs)){
	opt_S0 <- oa_path$satellites[which(oa_path$year==opt_start_year_bs[o])]
	opt_D0 <- oa_path$debris[which(oa_path$year==opt_start_year_bs[o])]
	opt_t0 <- which(oa_path$year==opt_start_year_bs[o]) - 1
	new_beginning <- opt_start_year_bs[o]
	message(paste0("\nGenerating optimal path beginning in ", new_beginning," for draw ",b,"...\n"))
	sink("log.solve.txt",append=TRUE)
	opt_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
	opt_tps_path <- tps_path_gen(opt_S0,opt_D0,opt_t0,R_start,R_start_year,R_frac,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
	sink()
	opt_path_list[[o]] <- cbind(year=seq(from=opt_start_year_bs[o],by=1,length.out=(T-opt_t0)),opt_tps_path)
}
opt_path <- rbindlist(opt_path_list)

	#############################################################################
	# 4. Draw plots, write output
	#############################################################################

	OA_OPT_full <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
	OA_OPT_full <- merge(OA_OPT_full,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"),all=TRUE)

	OA_OPT_full <- merge(OA_OPT_full,econ_series,by=c("year"),all=TRUE)

	selected_years <- intersect(which(OA_OPT_full$year>=start_year),which(OA_OPT_full$year<=end_year))
	OA_OPT <- OA_OPT_full[selected_years,]
	OA_OPT$bootstrap_draw <- b

	# Price of Anarchy in terms of collision risk. 1 represents no loss to anarchy, larger numbers show larger losses from anarchy.
	OA_OPT$riskPoA <- (OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt)*(OA_OPT$satellites.opt/OA_OPT$satellites.oa)
	# Price of Anarchy in terms of flow welfare. 1 : no present gains or losses to anarchy, >1 : present losses to anarchy, <1 : present gains to anarchy.
	OA_OPT$flowWelfPoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa 
	# Price of Anarchy in terms of NPV of welfare. 1 : no permanent gains or losses to anarchy, >1 : permanent losses to anarchy, <1 : permanent gains to anarchy.
	OA_OPT$NPVPoA <- OA_OPT$fleet_vfn_path.opt/OA_OPT$fleet_vfn_path.oa 

	# Since we're using aggregate data we need to divide by the number of satellites to get things into per-satellite units. Dividing by open access # of satellites puts everything relative to the "initial condition of open access".
	OA_OPT$flow_welfare_loss <- (OA_OPT$fleet_flowv.oa/OA_OPT$satellites.oa - OA_OPT$fleet_flowv.opt/OA_OPT$satellites.opt)*norm_const
	OA_OPT$npv_oa_welfare <- (OA_OPT$fleet_vfn_path.oa/OA_OPT$satellites.oa)*norm_const
	OA_OPT$npv_opt_welfare <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$satellites.opt)*norm_const
	OA_OPT$npv_welfare_loss <- (OA_OPT$npv_oa_welfare - OA_OPT$npv_opt_welfare)
	OA_OPT$npv_welfare_gain <- (OA_OPT$npv_opt_welfare - OA_OPT$npv_oa_welfare)

	F_over_horizon <- F[1:nrow(OA_OPT)]
	OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa/OA_OPT$satellites.oa - OA_OPT$collision_rate.opt/OA_OPT$satellites.opt)*F_over_horizon*norm_const*1e+9/OA_OPT$satellites.oa # 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.
	fe_eqm_next_ts <- data.frame(year=unique(OA_OPT$year),fe_eqm_next=fe_eqm[2:(length(unique(OA_OPT$year))+1)],launch_con=launch_constraint[2:(length(unique(OA_OPT$year))+1)])
	OA_OPT <- merge(OA_OPT,fe_eqm_next_ts)
	OA_OPT$num_destr_asat[is.na(OA_OPT$num_destr_asat)] <- 0

	OA_OPT$one_period_launch_deviation <- rep(-1,length=nrow(OA_OPT))
	for(i in 1:nrow(OA_OPT)){
		OA_OPT$one_period_launch_deviation[i] <- oa_deviation(OA_OPT$satellites.opt[i], OA_OPT$debris.opt[i], OA_OPT$fe_eqm_next[i], OA_OPT$launch_con[i], OA_OPT$num_destr_asat[i])
	}

	OA_OPT$one_period_sat_deviation <- S_(OA_OPT$one_period_launch_deviation,OA_OPT$satellites.opt,OA_OPT$debris.opt)
	OA_OPT$one_period_deb_deviation <- D_(OA_OPT$one_period_launch_deviation,OA_OPT$satellites.opt,OA_OPT$debris.opt,OA_OPT$num_destr_asat)

	shift_up <- function(x,optx) {
		output <- c(optx[1],x[1:(length(x)-1)])
		return(output)
	}

	deviation_dfrm <- data.frame(year=OA_OPT$year, start_time.opt=OA_OPT$start_time.opt,
						sats_dev=OA_OPT$one_period_sat_deviation, optS=OA_OPT$satellites.opt,
						debs_dev=OA_OPT$one_period_deb_deviation, optD=OA_OPT$debris.opt, 
						launch_dev=OA_OPT$one_period_launch_deviation, opt_launch=OA_OPT$launches.opt,
						optL=OA_OPT$collision_rate.opt/OA_OPT$satellites.opt, 
						F_over_horizon=OA_OPT$costs.opt)

	deviation_dfrm <- ddply(deviation_dfrm, ~start_time.opt, transform, sats = shift_up(sats_dev,optS), debs = shift_up(debs_dev,optD))
	deviation_dfrm <- ddply(deviation_dfrm, ~start_time.opt, transform, L_dev = L(sats_dev,debs_dev)/sats_dev)
	deviation_dfrm <- ddply(deviation_dfrm, ~start_time.opt, transform, excess_L = L_dev-optL)

	deviation_dfrm$excess_L[which(deviation_dfrm$excess_L<1e-16)] <- 0

	deviation_dfrm$opt_dev_tax_path <- (deviation_dfrm$excess_L*deviation_dfrm$F_over_horizon*norm_const*1e+9)/deviation_dfrm$optS
	# this is the tax that would deter a one-period open access deviation from a given path

	# 1e+9 scales to units of dollars from units of billion dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

	OA_OPT <- merge(OA_OPT,deviation_dfrm[,c("year","start_time.opt","excess_L","opt_dev_tax_path")],by=c("year","start_time.opt"))

	total_loop_wall_time <- round((proc.time()[3] - total_time)/60,3)
	opt_models_time_taken <- round((opt_end_time - opt_start_time)/60,3)
	opt_models_time_proportion <- round(opt_models_time_taken/total_loop_wall_time,3)*100
	message(paste0("\n Done. Total loop wall time: ",total_loop_wall_time," minutes. Proportion spent on calculating optimal models: ",opt_models_time_proportion,"%.\n"))
	if(b==1) {
		OA_OPT_bootstrap <- OA_OPT
	}
	if(b>1){
		OA_OPT_bootstrap <- rbind(OA_OPT_bootstrap, OA_OPT)
	}
}

write.csv(OA_OPT_bootstrap,file="../data/bootstrapped_simulation.csv")
