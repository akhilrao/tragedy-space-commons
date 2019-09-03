#####
# 1. This block is identical to the first block of removal_comparison_figures.r
# TODO: Replace this with block with a single script or functional block to achieve the same goal. Since similar things are being done here, in removal_comparison_figures, and main_model_figures, it should be doable.
#####

OA_OPT_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
OA_OPT_no_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))

#####

keep_cols <- c("year","launches.oa","launches.opt","satellites.oa","satellites.opt","debris.oa","debris.opt","fleet_vfn_path.oa","fleet_vfn_path.opt","collision_rate.oa","collision_rate.opt","start_time.opt","R_frac.opt","payloads_in_orbit","launch_successes","debris","risk.x","costs.opt")

OA_OPT_removal <- OA_OPT_removal[,keep_cols]
OA_OPT_no_removal <- OA_OPT_no_removal[,keep_cols]

OA_OPT <- merge(OA_OPT_removal,OA_OPT_no_removal,by=c("year","launch_successes","debris","risk.x","start_time.opt","payloads_in_orbit","costs.opt"),suffix=c(".rem",".norem"))

OA_OPT$riskPoA.rem <- (OA_OPT$collision_rate.oa.rem/OA_OPT$collision_rate.opt.rem)*(OA_OPT$satellites.opt.rem/OA_OPT$satellites.oa.rem)
OA_OPT$riskPoA.norem <- (OA_OPT$collision_rate.oa.norem/OA_OPT$collision_rate.opt.norem)*(OA_OPT$satellites.opt.norem/OA_OPT$satellites.oa.norem)
# Price of Anarchy in terms of NPV of welfare. 1 : no permanent gains or losses to anarchy, >1 : permanent losses to anarchy, <1 : permanent gains to anarchy.
OA_OPT$NPVPoA.rem <- OA_OPT$fleet_vfn_path.opt.rem/OA_OPT$fleet_vfn_path.oa.rem
OA_OPT$NPVPoA.norem <- OA_OPT$fleet_vfn_path.opt.norem/OA_OPT$fleet_vfn_path.oa.norem
OA_OPT$NPVPoA.oaremvnorem <- OA_OPT$fleet_vfn_path.oa.rem/OA_OPT$fleet_vfn_path.oa.norem
OA_OPT$NPVPoA.optremvnorem <- OA_OPT$fleet_vfn_path.opt.rem/OA_OPT$fleet_vfn_path.opt.norem
OA_OPT$NPVPoA.oaremvoptnorem <- OA_OPT$fleet_vfn_path.oa.rem/OA_OPT$fleet_vfn_path.opt.norem

# Since we're using aggregate data we need to divide by the number of satellites to get things into per-satellite units. To be consistent across scenarios, we use "satellites.oa.norem" as the division factor everywhere.
OA_OPT$npv_oa_welfare.rem <- (OA_OPT$fleet_vfn_path.oa.rem/OA_OPT$satellites.oa.norem)*norm_const
OA_OPT$npv_opt_welfare.rem <- (OA_OPT$fleet_vfn_path.opt.rem/OA_OPT$satellites.oa.norem)*norm_const
OA_OPT$npv_welfare_loss.rem <- (OA_OPT$npv_oa_welfare.rem - OA_OPT$npv_opt_welfare.rem)
OA_OPT$npv_welfare_gain.rem <- (OA_OPT$npv_opt_welfare.rem - OA_OPT$npv_oa_welfare.rem)

OA_OPT$npv_oa_welfare.norem <- (OA_OPT$fleet_vfn_path.oa.norem/OA_OPT$satellites.oa.norem)*norm_const
OA_OPT$npv_opt_welfare.norem <- (OA_OPT$fleet_vfn_path.opt.norem/OA_OPT$satellites.oa.norem)*norm_const
OA_OPT$npv_welfare_loss.norem <- (OA_OPT$npv_oa_welfare.norem - OA_OPT$npv_opt_welfare.norem)
OA_OPT$npv_welfare_gain.norem <- (OA_OPT$npv_opt_welfare.norem - OA_OPT$npv_oa_welfare.norem)

OA_OPT$npv_oa_welfare_gain_from_rem <- (OA_OPT$fleet_vfn_path.oa.rem - OA_OPT$fleet_vfn_path.oa.norem)

F_over_horizon <- OA_OPT$costs.opt

OA_OPT$opt_tax_path.rem <- (OA_OPT$collision_rate.oa.rem/OA_OPT$satellites.oa.rem - OA_OPT$collision_rate.opt.rem/OA_OPT$satellites.opt.rem)*F_over_horizon*norm_const*1e+9/OA_OPT$satellites.oa.norem 
OA_OPT$opt_tax_path.norem <- (OA_OPT$collision_rate.oa.norem/OA_OPT$satellites.oa.norem - OA_OPT$collision_rate.opt.norem/OA_OPT$satellites.opt.norem)*F_over_horizon*norm_const*1e+9/OA_OPT$satellites.oa.norem 

# 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters (the rescaling makes the value function iteration better-behaved). We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

ss_rows <- which(OA_OPT$start_time.opt==-1)

OA_OPT$start_year <- rep(-1,length=nrow(OA_OPT))
OA_OPT_SS <- OA_OPT[ss_rows,c("year","npv_opt_welfare.rem","npv_opt_welfare.norem")]
colnames(OA_OPT_SS)[2:3] <- c("ss_npv_opt_welfare.rem","ss_npv_opt_welfare.norem")
OA_OPT <- OA_OPT[-ss_rows,]

# relabel start_time.opt values from integer labels to the appropriate years
for(s in 1:length(unique(OA_OPT$start_time.opt))){
	OA_OPT$start_year[OA_OPT$start_time.opt==sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]] <- unique(OA_OPT$year)[sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]+1]
}

OA_OPT <- merge(OA_OPT,OA_OPT_SS,by=c("year"))

#####
# 2. This block also copies code from removal_comparison_figures.r
#####

boa_base_dfrm <- OA_OPT[which(OA_OPT$year>=2030),c("npv_welfare_gain.rem","npv_welfare_gain.norem","start_year","year")]

coi_base_dfrm <- boa_base_dfrm[which(boa_base_dfrm$start_year>2010),]

# npv_welfare_loss.rem is the cost of inaction (forgone gain from switching to opt mgmt in t instead of 2020, where t > 2020) with removal, .norem is without removal
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, npv_welfare_loss.norem=(npv_welfare_gain.norem[which(start_year==2020)]-npv_welfare_gain.norem), npv_welfare_loss.rem=(npv_welfare_gain.rem[which(start_year==2020)]-npv_welfare_gain.rem) )
# coi_effect_of_removal is the change in the cost of inaction due to introducing debris removal
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal=npv_welfare_loss.rem-npv_welfare_loss.norem)
# coi_effect_of_removal_pc is that change in cost of inaction expressed as a percentage of the gain from switching in 2020 with removal. "coi_effect_of_removal_pc=(coi_effect_of_removal/npv_welfare_gain.rem[which(start_year==2020)])*100" is the answer to, "supposing we had switched in 2020, with removal expected to come online in R_start_year, how much does the cost-of-inaction of switching in year t change when removal comes online?" so coi_effect_of_removal_pc = 7% says, "the introduction of debris removal in R_start_year increases the cost-of-inaction from switching in t by 7%, relative to having switched in 2020 with removal coming online in R_start_year"
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal_pc=(coi_effect_of_removal/npv_welfare_gain.rem[which(start_year==2020)])*100 )


# write out the files!
write.csv(coi_base_dfrm,file=paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_coi_base_dfrm.csv"))
