##### Function to calculate tax path given a dataframe. Meant to be called in main_model_figures.r, and for discount rate sensitivity analyisis. Takes a dataframe named OA_OPT as input, returns it as output with new columns added.

calc_tax_path <- function(OA_OPT,...) {
	# Price of Anarchy in terms of collision risk. 1 represents no loss to anarchy, larger numbers show larger losses from anarchy.
	OA_OPT$riskPoA <- (OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt)*(OA_OPT$satellites.opt/OA_OPT$satellites.oa)
	# Price of Anarchy in terms of flow welfare. 1 : no present gains or losses to anarchy, >1 : present losses to anarchy, <1 : present gains to anarchy.
	OA_OPT$flowWelfPoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa 
	# Price of Anarchy in terms of NPV of welfare from the fleet and then per-satellite. 1 : no permanent gains or losses to anarchy, >1 : permanent losses to anarchy, <1 : permanent gains to anarchy.
	OA_OPT$NPVPoA <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$fleet_vfn_path.oa)
	OA_OPT$NPVPoA_sat <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$fleet_vfn_path.oa)*(OA_OPT$satellites.oa/OA_OPT$satellites.opt)

	# Since we're using aggregate data we need to divide by the number of satellites to get the dollar values into per-fleet units. Otherwise, the dollar values are scaled by 2x the number of satellites rather than 1x -- 1x from the form of the pre-value function used in computation, and 1x from the aggregate dollar amounts used for calibration. These dollar amounts are in units of billion USD.
	OA_OPT$flow_welfare_loss <- (OA_OPT$fleet_flowv.oa/OA_OPT$satellites.oa - OA_OPT$fleet_flowv.opt/OA_OPT$satellites.opt)*norm_const
	#OA_OPT$npv_oa_welfare <- (OA_OPT$fleet_vfn_path.oa/OA_OPT$satellites.oa)*norm_const
	#OA_OPT$npv_opt_welfare <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$satellites.opt)*norm_const
	OA_OPT$npv_oa_welfare <- (OA_OPT$fleet_vfn_path.oa/OA_OPT$satellites.oa)*norm_const
	OA_OPT$npv_opt_welfare <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$satellites.oa)*norm_const
	OA_OPT$npv_welfare_loss <- (OA_OPT$npv_oa_welfare - OA_OPT$npv_opt_welfare)
	OA_OPT$npv_welfare_gain <- (OA_OPT$npv_opt_welfare - OA_OPT$npv_oa_welfare)

	# A tax which prevents returning to BAU from optimal path
	F_over_horizon <- OA_OPT$costs.opt
	#OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa/OA_OPT$satellites.oa - OA_OPT$collision_rate.opt/OA_OPT$satellites.opt)*(F_over_horizon*norm_const/OA_OPT$satellites.opt)*1e+9 
	OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa/OA_OPT$satellites.oa - OA_OPT$collision_rate.opt/OA_OPT$satellites.opt)*(F_over_horizon*norm_const/OA_OPT$satellites.oa)*1e+9 

	OA_OPT$num_destr_asat[is.na(OA_OPT$num_destr_asat)] <- 0

	ss_rows <- which(OA_OPT$start_time.opt==-1)

	OA_OPT$start_year <- rep(-1,length=nrow(OA_OPT))
	OA_OPT_SS <- OA_OPT[ss_rows,c("year","npv_opt_welfare")]
	colnames(OA_OPT_SS)[2] <- "ss_npv_opt_welfare"
	OA_OPT <- OA_OPT[-ss_rows,]

	# Loop to convert start_time.opt codes into start_year with year-labels
	for(s in 1:length(unique(OA_OPT$start_time.opt))){
		OA_OPT$start_year[OA_OPT$start_time.opt==sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]] <- unique(OA_OPT$year)[sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]+1]
	}

	OA_OPT <- merge(OA_OPT,OA_OPT_SS,by=c("year"))

	# A tax which prevents OA jump off of an optimal path
	fe_eqm_next_ts <- data.frame(year=unique(OA_OPT$year),fe_eqm_next=fe_eqm[2:(length(unique(OA_OPT$year))+1)],launch_con=launch_constraint[2:(length(unique(OA_OPT$year))+1)])
	OA_OPT <- merge(OA_OPT,fe_eqm_next_ts)

	OA_OPT$one_period_launch_deviation <- rep(-1,length=nrow(OA_OPT))
	for(i in 1:nrow(OA_OPT)){
		OA_OPT$one_period_launch_deviation[i] <- oa_deviation(OA_OPT$satellites.opt[i], OA_OPT$debris.opt[i], OA_OPT$fe_eqm_next[i], OA_OPT$launch_con[i], OA_OPT$num_destr_asat[i])
	}

	OA_OPT$one_period_sat_deviation <- S_(OA_OPT$one_period_launch_deviation,OA_OPT$satellites.opt,OA_OPT$debris.opt)
	OA_OPT$one_period_deb_deviation <- D_(OA_OPT$one_period_launch_deviation,OA_OPT$satellites.opt,OA_OPT$debris.opt,OA_OPT$num_destr_asat)

	deviation_dfrm <- data.frame(year=OA_OPT$year, start_year=OA_OPT$start_year,
						sats_orig=OA_OPT$one_period_sat_deviation, optS=OA_OPT$satellites.opt,
						debs_orig=OA_OPT$one_period_deb_deviation, optD=OA_OPT$debris.opt, 
						launch_dev=OA_OPT$one_period_launch_deviation, opt_launch=OA_OPT$launches.opt,
						optL=OA_OPT$collision_rate.opt/OA_OPT$satellites.opt, 
						F_over_horizon=OA_OPT$costs.opt)

	OA_OPT$collision_rate.opt[which(deviation_dfrm$start_year==2020)]/OA_OPT$satellites.opt[which(deviation_dfrm$start_year==2020)]

	deviation_dfrm <- ddply(deviation_dfrm, ~start_year, transform, sats = shift_up(sats_orig,optS), debs = shift_up(debs_orig,optD))
	deviation_dfrm <- ddply(deviation_dfrm, ~start_year, transform, L_dev = L(sats,debs)/sats)
	deviation_dfrm <- ddply(deviation_dfrm, ~start_year, transform, excess_L = L_dev-optL)

	deviation_dfrm$excess_L[which(deviation_dfrm$excess_L<1e-16)] <- 0

	deviation_dfrm$opt_dev_tax_path <- (deviation_dfrm$excess_L*deviation_dfrm$F_over_horizon*norm_const*1e+9)/deviation_dfrm$optS
	# this is the tax that would deter a one-period open access deviation from a given path

	# 1e+9 scales to units of dollars from units of billion dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

	OA_OPT <- merge(OA_OPT,deviation_dfrm[,c("year","start_year","excess_L","opt_dev_tax_path")],by=c("year","start_year"))

	OA_OPT_full <- OA_OPT
	#OA_OPT <- OA_OPT[-which(OA_OPT_full$year>2041),] # this line truncates back to a specific period, seems to have caused some issues in the discount rate counterfactual... doesn't seem necessary?
	return(OA_OPT)
}

shift_up <- function(x,optx) {
	output <- c(optx[1],x[1:(length(x)-1)])
	return(output)
}