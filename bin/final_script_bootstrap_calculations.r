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

cat(paste0("\n\nBeginning bootstrap draw ",b))
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
	S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay #guess for BW parameterization
	V_T <- p[T]*S_T
	vguess <- matrix(V_T,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
	lpguess <- matrix(0,nrow=S_gridsize_opt,ncol=D_gridsize_opt)
	gridpanel <- grid_to_panel(opt_gridlist,lpguess,vguess)

	# initialize solver output list
	opt_dvs_output <- list()

	# run path solver
cat(paste0("\nCalculating optimal models for draw ",b,"..."))
sink("log.solve.txt", append=TRUE)
	opt_dvs_output <- suppressWarnings(opt_pvfn_path_solver(opt_dvs_output,gridpanel,S_gridsize_opt,D_gridsize_opt,opt_gridlist,asats,T,p,F,ncores=ncores))
sink()

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
	S_T <- (S_T_1 - L(S_T_1,D_T_1))*avg_sat_decay #guess for BW parameterization
	V_T <- p[T]*S_T
	vguess <- matrix(V_T,nrow=gridsize,ncol=gridsize)
	lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
	gridpanel <- grid_to_panel(oa_gridlist,lpguess,vguess)

	# initialize solver output list
	oa_dvs_output <- list()

	# run path solver
cat(paste0("\nCalculating open access models for draw ",b,"..."))
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
cat(paste0("\nGenerating open access path for draw ",b,"..."))
sink("log.solve.txt", append=TRUE)
	oa_tps_path <- tps_path_gen(S0,D0,0,R_start,R_start_year,R_frac,p,F,oa_pvfn_path,asats,launch_constraint,oa_grid_lookup,ncores=ncores,OPT=0,linear_policy_interp=0)
sink()
	oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(oa_tps_path)),oa_tps_path)

	### optimal paths
	opt_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
cat(paste0("\nGenerating optimal path for draw ",b,"..."))
sink("log.solve.txt", append=TRUE)
	opt_tps_path <- tps_path_gen(S0,D0,0,R_start,R_start_year,R_frac,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
sink()
	opt_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(opt_tps_path)),opt_tps_path)

	#############################################################################
	# 4. Draw plots, write output
	#############################################################################

	OA_OPT_full <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
	OA_OPT_full <- merge(OA_OPT_full,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"),all=TRUE)

	OA_OPT_full <- merge(OA_OPT_full,econ_series,by=c("year"),all=TRUE)

	selected_years <- intersect(which(OA_OPT_full$year>=start_year),which(OA_OPT_full$year<=end_year))
	OA_OPT <- OA_OPT_full[selected_years,]
	OA_OPT$bootstrap_draw <- b

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

	cat(paste0("\n Done. Total loop wall time: ",round(proc.time()[3] - total_time,3)/60," minutes"))
	if(b==1) {
		OA_OPT_bootstrap <- OA_OPT
	}
	if(b>1){
		OA_OPT_bootstrap <- rbind(OA_OPT_bootstrap, OA_OPT)
	}
}

write.csv(OA_OPT_bootstrap,file="../data/bootstrapped_simulation.csv")
