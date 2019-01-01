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
source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

#############################################################################
# Calibration
#############################################################################

risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

T <- nrow(econ_series)
upper <- 1e15
#upper_seq <- seq(from=1,by=1,length=T)

start_year <- 2006
S0 <- observed_time_series$payloads_in_orbit[which(observed_time_series$year==start_year)]
D0 <- observed_time_series$debris[which(observed_time_series$year==start_year)]
avg_sat_decay <- mean(observed_time_series$payloads_decayed[which(observed_time_series$year>=start_year)]/observed_time_series$payloads_in_orbit[which(observed_time_series$year>=start_year)])

aS <- risk_cal[1,2]
aD <- risk_cal[2,2]
aSS <- risk_cal[3,2]
aSD <- risk_cal[4,2]
aDD <- risk_cal[5,2]

aDDbDD <- deblom_cal[7,2]
bSS <- deblom_cal[5,2]
bSD <- deblom_cal[6,2]
d <- deblom_cal[2,2]
Z_coef <- deblom_cal[1,2]
m <- deblom_cal[3,2]
asat_coef <- deblom_cal[4,2]
asats <- observed_time_series$num_destr_asat[which(observed_time_series$year>=start_year)]

discount_rate <- 0.05#0.6660755
discount_fac <- 1/(1+discount_rate)

fe_eqm <- econ_coefs[1,2] + econ_coefs[2,2]*econ_series$r_s + econ_coefs[3,2]*econ_series$Ft_Ft
fe_eqm <- observed_time_series$risk[which(observed_time_series$year>2005)]
#p <- rep(1,length=length(fe_eqm))
# function to calculate sequence of economic launch costs recursively
F_calc <- function(a_1,a_2,a_3,F_1,pi_t,risk,...) {
	F <- rep(0,length=length(risk))
	F[1] <- F_1
	for(i in 1:(length(risk)-1)) {
		F[i+1] <- (a_2*pi_t[i] + a_3*F[i])/(risk[i+1] - a_1)
	}
	return(F)
}

#F <- F_calc(a_1=econ_coefs[1,2],a_2=econ_coefs[2,2],a_3=econ_coefs[3,2],F_1=(econ_series$Ft_Ft[1]*econ_series$tot_cost[1]),pi_t=econ_series$tot_rev,risk=fe_eqm)

physecon <- merge(observed_time_series,econ_series,by=c("year","risk"))
p_per_S <- physecon$tot_rev/physecon$payloads_in_orbit
F_per_S <- physecon$tot_cost/physecon$payloads_in_orbit
# fe_eqm <- p/F - discount_rate*c(F[1],F[-length(F)])/F

observed_launches <- observed_time_series$launch_successes[which(observed_time_series$year>=start_year)]+observed_time_series$launch_failures[which(observed_time_series$year>=start_year)]
launch_constraint <- 1000*cummax(observed_launches)
#launch_constraint[25:length(launch_constraint)] <- 0.05*max(observed_launches)

#############################################################################
# Open access time path
#############################################################################

oa_path <- oa_tsgen(S0,D0,T,fe_eqm,launch_constraint,asats)

oa_path$time <- start_year+oa_path$time
colnames(oa_path)[which(colnames(oa_path)=="time")] <- "year"

##### compare fit of generated series against actual
oa_merged <- merge(oa_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

oa_merged <- cbind(oa_merged,fe_eqm=fe_eqm)

oa_fitcomp_base <- ggplot(data=oa_merged,aes(x=year))
oa_launch_fit <- oa_fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
oa_sat_fit <- oa_fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
oa_deb_fit <- oa_fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
oa_risk_fit <- oa_fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

dev.new()
#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(oa_launch_fit,oa_sat_fit,oa_risk_fit,oa_deb_fit,ncol=2)
#dev.off()
View(oa_merged)
#############################################################################
# Optimal time path
#############################################################################

# build grid, generate guesses, initialize dynamic_vfi_solver output list
gridsize <- 32
gridlist <- build_grid(gridmin=0, gridmax=35000, gridsize, cheby=1)
### shift grid down to zero if chebysheving it moved it up - should be unnecessary with expanded chebysehv array (secant factor)
if(min(gridlist$base_piece)>0) {
	gridlist$base_piece <- gridlist$base_piece - min(gridlist$base_piece)
	gridlist$igrid <- gridlist$igrid - min(gridlist$igrid)
}
vguess <- matrix(gridlist$igrid$sats,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
dvs_output <- list()

p <- p_per_S
F <- F_per_S
# Check that path solver works in all periods
dvs_output <- policy_function_path_solver(gridpanel,gridlist,asats,T,p,F,ncores=4)

# bind the list of solved policies into a long dataframe
policy_path <- rbindlist(dvs_output)

# compare tps time series to solved time series
grid_lookup <- data.frame(sats=policy_path$satellites,debs=policy_path$debris,F=policy_path$F)

tps_path <- tps_opt_path(S0,D0,p,F,policy_path,asats,launch_constraint,grid_lookup,launchcon_seq=launch_constraint,ncores=4)

opt_path <- cbind(year=seq(from=start_year,by=1,length.out=nrow(tps_path)),tps_path)


##### compare fit of generated series against actual
opt_merged <- merge(opt_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

opt_merged <- cbind(opt_merged,fe_eqm=fe_eqm)

opt_fitcomp_base <- ggplot(data=opt_merged,aes(x=year))
opt_launch_fit <- opt_fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
opt_sat_fit <- opt_fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
opt_deb_fit <- opt_fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
opt_risk_fit <- opt_fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

dev.new()
#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(opt_launch_fit,opt_sat_fit,opt_risk_fit,opt_deb_fit,ncol=2)
#dev.off()
View(opt_merged)
#############################################################################
# Compare OA and OPT outcomes
#############################################################################

OA_OPT <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT <- merge(OA_OPT,econ_series,by=c("year"))
OA_OPT$oa_value <- p_per_S*OA_OPT$satellites.oa + (1-OA_OPT$collision_rate.oa)*F_per_S*OA_OPT$satellites.oa - F_per_S*OA_OPT$launches.oa
OA_OPT$fp_value <- p_per_S*(OA_OPT$satellites.oa) + (1-OA_OPT$collision_rate.oa)*F_per_S*(OA_OPT$satellites.oa) - F_per_S*(OA_OPT$launches.oa)
#OA_OPT$fp_value <- p_per_S*OA_OPT$satellites.opt + (1-OA_OPT$collision_rate.opt)*F_per_S*OA_OPT$satellites.opt


OA_OPT$riskPoA <- OA_OPT$collision_rate.opt/OA_OPT$collision_rate.oa
OA_OPT$PoA <- OA_OPT$fp_value/OA_OPT$oa_value
OA_OPT$welfare_loss <- OA_OPT$fp_value - OA_OPT$oa_value

oaoptcomp_base <- ggplot(data=OA_OPT,aes(x=year))

risk_comps <- oaoptcomp_base + geom_line(aes(y=riskPoA),size=0.85) +
							geom_hline(yintercept=1,linetype="dashed",color="blue") +
							ylab("Price of Anarchy for launch decisions") + xlab("year") + theme_minimal()

welf_loss <- oaoptcomp_base + geom_line(aes(y=welfare_loss),size=0.85) +
							ylab("Welfare loss from open access ($100m/year)") + xlab("year") + theme_minimal()
dev.new()
grid.arrange(welf_loss,risk_comps,ncol=2)
