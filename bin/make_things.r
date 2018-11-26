rm(list=ls())

library(rootSolve)
library(gridExtra)
library(ggplot2)
library(viridis)
source("equations.r")
source("simulation_functions.r")
source("simulation_algorithms.r")

#############################################################################
# Function definitions
#############################################################################



#############################################################################
# Calibration
#############################################################################

risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

# append fake steady state
# nfakes <- 100
# fakedata <- matrix(0,nrow=nfakes,ncol=dim(observed_time_series)[2])
# lastrow <- as.numeric(observed_time_series[dim(observed_time_series)[1],])
# for(i in 1:nfakes) {
# 	fakedata[i,] <- lastrow
# 	fakedata[i,2] <- 2015+i
# 	fakedata[i,ncol(fakedata)] <- ifelse(i<nfakes/4,lastrow[ncol(fakedata)] + 0.0025*i,max(fakedata[,ncol(fakedata)]))
# }
# colnames(fakedata) <- colnames(observed_time_series)
# observed_time_series <- rbind(observed_time_series,fakedata)
# T <- nfakes+11

T <- nrow(econ_series)
upper <- 1e15
upper_seq <- seq(from=1,by=1,length=T)

start_year <- 2005
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
launch_constraint <- cummax(observed_launches)
#launch_constraint[25:length(launch_constraint)] <- 0.05*max(observed_launches)

#############################################################################
# Open access time path
#############################################################################

oa_path <- oa_tsgen(S0,D0,T,fe_eqm,launch_constraint,asats)

oa_path$time <- start_year+oa_path$time
colnames(oa_path)[which(colnames(oa_path)=="time")] <- "year"

##### compare fit of generated series against actual
oa_merged <- merge(oa_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

oa_merged <- cbind(oa_merged,fe_eqm=c(fe_eqm[1],fe_eqm))

fitcomp_base <- ggplot(data=oa_merged,aes(x=year))
launch_fit <- fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
sat_fit <- fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
deb_fit <- fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
risk_fit <- fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(launch_fit,sat_fit,risk_fit,deb_fit,ncol=2)
#dev.off()

#############################################################################
# Optimal time path
#############################################################################

### function to begin an optimal launch sequence at a given time
simulate_optimal_path <- function(start_year,S0,D0,pi_t,F_t,discount_rate,launch_con,asats,...) {
	T <- 500
	fe_eqm <- pi_t/F_t - discount_rate
	fe_eqm_inf <- c(fe_eqm,rep(fe_eqm[length(fe_eqm)],length=(T-length(fe_eqm))))
	asats_inf <- c(asats,rep(0,length=(T-length(fe_eqm)-1)))
	launch_constraint_inf <- c(launch_constraint,rep(launch_constraint[length(launch_constraint)],length=(T-length(fe_eqm)-1)))

	opt_path <- fp_tsgen(S0,D0,T,fe_eqm_inf,launch_constraint_inf,asats_inf)

	check <- seriesgen_ts(rep(1,length=T),S0,D0,T,asats_inf)

	opt_path$time <- start_year+opt_path$time
	colnames(opt_path)[which(colnames(opt_path)=="time")] <- "year"
	return(opt_path)
}

p <- p_per_S
F <- F_per_S
opt_path <- simulate_optimal_path(start_year,S0,D0,econ_series$tot_rev,econ_series$tot_cost,discount_rate,launch_constraint,asats)



##### compare fit of generated series against actual
opt_merged <- merge(opt_path,observed_time_series,by=c("year"),suffixes=c(".sim",".obs"))

opt_merged <- cbind(opt_merged,fe_eqm=c(fe_eqm[1],fe_eqm))

fitcomp_base <- ggplot(data=opt_merged,aes(x=year))
launch_fit <- fitcomp_base + geom_line(aes(y=launches),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=launch_successes),size=1) +
							ylab("Satellites launched") + theme_minimal() +
							ggtitle("Simulated series (blue) vs. observed (black)")
sat_fit <- fitcomp_base + geom_line(aes(y=satellites),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=payloads_in_orbit),size=1) +
							ylab("Satellites in LEO") + theme_minimal() +
							ggtitle("")
deb_fit <- fitcomp_base + geom_line(aes(y=debris.sim),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=debris.obs),size=1) +
							ylab("Debris in LEO") + xlab("year") + theme_minimal()
risk_fit <- fitcomp_base + geom_line(aes(y=collision_rate),linetype="dashed",color="blue",size=0.85) +
							geom_line(aes(y=risk),size=1) +
							geom_line(aes(y=fe_eqm),size=1,linetype="dotted") +
							ylab("Collision risk in LEO") + xlab("year") + theme_minimal()

#png(width=700,height=700,filename="../images/open_access_historical_fit.png")
grid.arrange(launch_fit,sat_fit,risk_fit,deb_fit,ncol=2)
#dev.off()

#############################################################################
# Generate optimal policy functions along time path
#############################################################################


gridmin <- 0
gridmax <- 10000
ss_vguess <- matrix(0,nrow=25,ncol=25)
ss_lpguess <- matrix(0,nrow=25,ncol=25)

ss_vfi_solver(vguess,launch_pguess,gridmin,gridmax)



#############################################################################
# Compare OA and OPT outcomes
#############################################################################

OA_OPT <- merge(oa_path,opt_path,by=c("year"),suffixes=c(".oa",".opt"))
OA_OPT <- merge(OA_OPT,econ_series,by=c("year"))
OA_OPT$oa_value <- p_per_S*OA_OPT$satellites.oa + (1-OA_OPT$collision_rate.oa)*F_per_S*OA_OPT$satellites.oa - F_per_S*OA_OPT$launches.oa
OA_OPT$fp_value <- p_per_S*(OA_OPT$satellites.oa+1) + (1-OA_OPT$collision_rate.oa)*F_per_S*(OA_OPT$satellites.oa+1) - F_per_S*(OA_OPT$launches.oa+1)
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
grid.arrange(welf_loss,risk_comps,ncol=2)
