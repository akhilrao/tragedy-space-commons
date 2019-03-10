##### Script to set calibrated economic and physical parameters based on observed data and theory

# Read in data and parameter estimates
## bootstrapped parameters
risk_cal_set <- read.csv("../data/bootstrapped_risk_eqn_coefs.csv")
test <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal_set <- read.csv("../data/bootstrapped_debris_lom_coefs.csv")
## non-bootstrapped values
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
implied_econ_series <- read.csv("../data/implied_costs.csv")
satlom_cal <- read.csv("../data/calibrated_satellite_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")
MS_proj_rev <- read.csv("../data/avg_econ_return.csv")
MS_proj_total <- read.csv("../data/avg_econ_total.csv")

# restrict risk_cal parameters to values where SD > 0 (both collision risk couplings active)
accepted_risk_cal_set <- risk_cal_set[which(risk_cal_set$SD>0),]

# the physical parameters are draws from the bootstrap world's conditional distribution, parameters(risk) and parameters(debris|risk)
set.seed(501)
start_loc <- sample(c(1:(nrow(risk_cal_set)-B)),size=1)
risk_cal_set_B <- accepted_risk_cal_set[start_loc:(start_loc+B),-1]
deblom_cal_set_B <- deblom_cal_set[start_loc:(start_loc+B),-1]
# bootstrap_grid <- matrix(-1,nrow=B*B,ncol=sum(ncol(risk_cal_set_B),ncol(deblom_cal_set_B)))
# for(i in 1:B) {
# 	for(j in 1:B) {
# 		bootstrap_grid[(j+(i-1)*B),1] <- risk_cal_set_B[j,1]
# 		bootstrap_grid[(j+(i-1)*B),2] <- risk_cal_set_B[j,2]
# 		bootstrap_grid[(j+(i-1)*B),3] <- deblom_cal_set_B[i,1]
# 		bootstrap_grid[(j+(i-1)*B),4] <- deblom_cal_set_B[i,2]
# 		bootstrap_grid[(j+(i-1)*B),5] <- deblom_cal_set_B[i,3]
# 		bootstrap_grid[(j+(i-1)*B),6] <- deblom_cal_set_B[i,4]
# 		bootstrap_grid[(j+(i-1)*B),7] <- deblom_cal_set_B[i,5]
# 	}
# }

bootstrap_grid <- cbind(risk_cal_set_B,deblom_cal_set_B)
bootstrap_grid <- data.frame(bootstrap_grid)
colnames(bootstrap_grid) <- c(colnames(risk_cal_set_B),colnames(deblom_cal_set_B))

# select bootstrap parameters
risk_cal <- bootstrap_grid[b,1:2]#accepted_risk_cal_set[b,]
deblom_cal <- bootstrap_grid[b,3:7]#deblom_cal_set[b,]

# Extend Morgan Stanley revenue and total value projections an additional 5 years, to avoid any end-of-horizon effects for a forecast out to 2040 (e.g. numerical distortions in steady-state value functions). The idea is to "project" out to 2050 using the mean annual growth rate of the Morgan Stanley projections, then truncate back to 2040 to avoid any end-of-horizon effects.
projection_start <- MS_proj_rev$Year[nrow(MS_proj_rev)]+1
projection_end <- 2050
revenue_mean_CAGR <- mean((MS_proj_rev[-1,2]/MS_proj_rev[-nrow(MS_proj_rev),2] - 1))
revenue_growth <- cumprod(rep(1+revenue_mean_CAGR,length=length(c(projection_start:projection_end))))
revenue_projection <- data.frame(Year=c(projection_start:projection_end),Revenues=(MS_proj_rev[nrow(MS_proj_rev),2]*revenue_growth))
total_mean_CAGR <- mean((MS_proj_total[-1,2]/MS_proj_total[-nrow(MS_proj_total),2] - 1))
total_growth <- cumprod(rep(1+total_mean_CAGR,length=length(c(projection_start:projection_end))))
total_projection <- data.frame(Year=c(projection_start:projection_end),Total=(MS_proj_total[nrow(MS_proj_total),2]*total_growth))
MS_proj_rev <- rbind(MS_proj_rev,revenue_projection)
MS_proj_total<- rbind(MS_proj_total,total_projection)

# Set parameters
start_year <- 2006
end_year <- 2040
selected_years <- which(observed_time_series$year>=start_year)
S0 <- observed_time_series$payloads_in_orbit[which(observed_time_series$year==start_year)]
D0 <- observed_time_series$debris[which(observed_time_series$year==start_year)]

# estimated parameterization
satlom_cal_names <- as.character(satlom_cal[,1])
satlom_cal <- data.frame(parameters=t(c(satlom_cal[,2])))
colnames(satlom_cal) <- satlom_cal_names
avg_sat_decay <- satlom_cal$payloads_in_orbit # corresponds to just over 30 years on orbit: on average 5 year mission time + 25 year post-mission disposal compliance. Value estimated from statistical model for satellite law of motion. 
# decay coefficient is currently represented as survival rate rather than decay rate.

aSS <- risk_cal$S2
aSD <- risk_cal$SD

aDDbDD <- 0#deblom_cal$D2
bSS <- deblom_cal$SSfrags
bSD <- deblom_cal$SDfrags
d <- 1-deblom_cal$debris
#Z_coef <- deblom_cal[1,2]
m <- deblom_cal$launch_successes
asat_coef <- deblom_cal$num_destr_asat

discount_rate <- 0.05
discount_fac <- 1/(1+discount_rate)

# The implied series has the "true" economic costs and returns implied by the measurement error equation. The procedure brings in future economic information from the calibration. Combine this with the observed economic series.
comb_econ_series <- merge(econ_series,implied_econ_series,by=c("year")) 

# Merge the economic and physical data. The final period gets truncated because the theory-adjustment formula uses up one year.
#One option is to pad the final observation with the last value in the series. For p, this value is observed. For F, this value is just a copy of the newly-penultimate value. This carries some implications for the measurement errors in the final period.
physecon <- merge(observed_time_series,comb_econ_series,by=c("year","risk"))
p <- physecon$pi_t
F <- physecon$F_hat

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
fe_eqm <- econ_coefs[1,2] + econ_coefs[2,2]*econ_series$r_s + econ_coefs[3,2]*econ_series$Ft_Ft # calculate the path of the OA eqm condition from the calibrated regression

fe_eqm_proj <- econ_coefs[1,2] + econ_coefs[2,2]*(p[-length(p)]/F[-length(F)]) + econ_coefs[3,2]*(F[-1]/F[-length(F)]) # the final period, assume underlying econ parameters are constant

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
