##### Script to set calibrated economic and physical parameters based on observed data and theory

# Read in data and calibration estimates
risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
implied_econ_series <- read.csv("../data/implied_costs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")
MS_proj_rev <- read.csv("../data/avg_econ_return.csv")
MS_proj_total <- read.csv("../data/avg_econ_total.csv")

# Set parameters
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

discount_rate <- 0.05
discount_fac <- 1/(1+discount_rate)

# The implied series has the "true" economic costs and returns implied by the measurement error equation. The procedure brings in future economic information from the calibration. Combine this with the observed economic series.
comb_econ_series <- merge(econ_series,implied_econ_series,by=c("year")) 

# Merge the economic and physical data. The final period gets truncated because the theory-adjustment formula uses up one year.
#One option is to pad the final observation with the last value in the series. For p, this value is observed. For F, this value is just a copy of the newly-penultimate value. This carries some implications for the measurement errors in the final period.
physecon <- merge(observed_time_series,comb_econ_series,by=c("year","risk"))
p <- physecon$pi_t
#p <- c(p,econ_series$tot_rev[length(econ_series$tot_rev)])
F <- physecon$F_hat
#F <- c(F,F[length(F)])

# Append data with average case revenue projection from Morgan Stanley.
econ_proj <- merge(MS_proj_rev,MS_proj_total,by=c("Year"))
econ_proj$Costs <- econ_proj$Total - econ_proj$Revenues
econ_proj[,c(2,3,4)] <- econ_proj[,c(2,3,4)]/1000
p <- c(p,econ_proj$Revenues)
F <- c(F,econ_proj$Costs)

# Normalize pi_1 := 1. This helps convergence to a fixed tolerance because the numbers are smaller to begin with. In the end, we undo the normalization to get things back into the right units.
norm_const <- p[1]
p <- p/norm_const
F <- F/norm_const

##### Which method of generating the fe_eqm path is better? fe_eqm is used to generate the open access policy functions.

# Regression: This is what was estimated. It adjusts for measurement error.
# fe_eqm <- econ_coefs[1,2] + econ_coefs[2,2]*econ_series$r_s + econ_coefs[3,2]*econ_series$Ft_Ft # calculate the path of the OA eqm condition from the calibrated regression

# Observed: this is what was observed. it does not adjust for measurement error, but may fit the observed data better.
fe_eqm <- observed_time_series$risk[which(observed_time_series$year>2005)] # use the observed collision probability as the path of the OA eqm condition

fe_eqm_proj <- econ_coefs[1,2] + econ_coefs[2,2]*(p[-length(p)]/F[-length(F)]) + econ_coefs[3,2]*(F[-1]/F[-length(F)])

p <- p[-length(p)]
F <- F[-length(F)]

fe_eqm <- c(fe_eqm,fe_eqm_proj[-(1:length(fe_eqm))])

#####

# Set time-related values
T <- length(p)
asats <- observed_time_series$num_destr_asat[which(observed_time_series$year>=start_year)]
asats <- c(asats,rep(0,length=(T-length(asats))))

# Construct historical launch constraint. In each period, the number of possible launches are the largest number of launch attempts (failures+successes) observed until that period.
observed_launches <- observed_time_series$launch_successes[which(observed_time_series$year>=start_year)]+observed_time_series$launch_failures[which(observed_time_series$year>=start_year)]
obs_launch_constraint <- cummax(observed_launches)

# Relax historical launch constraint: have it increase at a linear trend rate
lc_data <- data.frame(period=seq(from=1,length.out=length(obs_launch_constraint)),limit=obs_launch_constraint)
lc_model <- lm(limit ~ period, data=lc_data)
lc_coef <- as.matrix(coef(lc_model))
lc_design <- as.matrix(data.frame(intercept=rep(1,length=(length(F)-length(obs_launch_constraint))),period=seq(from=11,length.out=(length(F)-length(obs_launch_constraint)))))
lc_proj <- lc_design%*%lc_coef

launch_constraint <- c(obs_launch_constraint,floor(lc_proj))

lc_dfrm <- data.frame(year=seq(from=start_year,length.out=length(launch_constraint)),observed_constraint=c(obs_launch_constraint, rep(NA,length=(length(launch_constraint)-length(obs_launch_constraint)))),projected_constraint=launch_constraint)

lc_path_base <- ggplot(data=lc_dfrm,aes(x=year))
lc_path <- lc_path_base + geom_line(aes(y=projected_constraint),linetype="dashed",color="blue",size=0.8) +
		geom_line(aes(y=observed_constraint),size=0.85) +							
		geom_vline(xintercept=2015,size=1,linetype="dashed") +
		ylab("Cumulative maximum launch attempts") + theme_minimal() +
		ggtitle("Observed and projected launch constraint")

lc_path

png(width=700,height=700,filename="../images/linear_trend_launch_constraint.png")
lc_path
dev.off()
