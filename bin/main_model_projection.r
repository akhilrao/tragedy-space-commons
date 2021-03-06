##### Script to generate main model projections for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Generate time paths from policy sequences
# 2. Write output

#############################################################################
# 1. Generate open access and optimal time paths
#############################################################################

path_compute_time <- proc.time()[3]

sink("log.project.txt",append=TRUE)

### open access paths
oa_grid_lookup <- data.frame(sats=oa_pvfn_path$satellites,debs=oa_pvfn_path$debris,F=oa_pvfn_path$F)
message("\nGenerating open access path...")
oa_tps_path <- tps_path_gen(S0,D0,0,R_start,R_start_year,R_frac,p,F,oa_pvfn_path,asats,launch_constraint,oa_grid_lookup,ncores=ncores,OPT=0,linear_policy_interp=0)
oa_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(oa_tps_path)),oa_tps_path)

### optimal paths
opt_path_list <- list()
for(o in 1:length(opt_start_year)){
	opt_S0 <- oa_path$satellites[which(oa_path$year==opt_start_year[o])]
	opt_D0 <- oa_path$debris[which(oa_path$year==opt_start_year[o])]
	opt_t0 <- which(oa_path$year==opt_start_year[o]) - 1
	new_beginning <- opt_start_year[o]
	message(paste0("\nGenerating optimal path beginning in ", new_beginning,"..."))
	opt_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
	opt_tps_path <- tps_path_gen(opt_S0,opt_D0,opt_t0,R_start,R_start_year,R_frac,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
	opt_path_list[[o]] <- cbind(year=seq(from=opt_start_year[o],by=1,length.out=(T-opt_t0)),opt_tps_path)
}
opt_path <- rbindlist(opt_path_list)

# opt_SS_S0 <- tail(opt_path$satellite[which(opt_path$start_time==0)])[6]
# opt_SS_D0 <- tail(opt_path$debris[which(opt_path$start_time==0)])[6]
# opt_SS_t0 <- 0
# message(paste0("\nGenerating optimal 'steady state' path beginning in 2006 at (S,D)=(", paste0(round(opt_SS_S0,digits=2),",",round(opt_SS_D0,digits=2)),")..."))
# opt_SS_grid_lookup <- data.frame(sats=opt_pvfn_path$satellites,debs=opt_pvfn_path$debris,F=opt_pvfn_path$F)
# opt_SS_tps_path <- tps_path_gen(opt_SS_S0,opt_SS_D0,opt_SS_t0,R_start,R_start_year,R_frac,p,F,opt_pvfn_path,asats,launch_constraint,opt_grid_lookup,ncores=ncores,OPT=1,linear_policy_interp=0)
# opt_SS_path <- cbind(year=seq(from=start_year,by=1,length.out=(T-opt_SS_t0)),opt_SS_tps_path)
# opt_SS_path$start_time <- -1

# opt_path <- rbind(opt_path,opt_SS_path)

sink()

path_compute_time <- round(proc.time()[3]-path_compute_time,3)
message("Time to generate paths: ",path_compute_time," seconds")


#############################################################################
# 2. Write output
#############################################################################

OA_OPT_full <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT_full <- merge(OA_OPT_full,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"),all=TRUE)

OA_OPT_full <- merge(OA_OPT_full,econ_series,by=c("year"),all=TRUE)

selected_years <- intersect(which(OA_OPT_full$year>=start_year),which(OA_OPT_full$year<=end_year))
OA_OPT <- OA_OPT_full[selected_years,]

if(counterfactual=="none") {
	write.csv(OA_OPT, file=paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
}
if(counterfactual=="avoidance") {
	write.csv(OA_OPT, file=paste0("../data/counterfactuals/collision_avoidance/",opt_start_year[1],"_cf_avoidance_aSS_",round(log(aSS),1),"_aSD_",round(log(aSD),1),"_simulation.csv"))
} 
if(counterfactual=="discount") {
	write.csv(OA_OPT, file=paste0("../data/counterfactuals/discount_rate/",opt_start_year[1],"_cf_discount_r_",discount_rate,"_simulation.csv"))
} 

setwd("../bin")
