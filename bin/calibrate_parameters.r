##### Script to set calibrated economic and physical parameters based on observed data and theory

# Read in data and calibration estimates
risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
implied_econ_series <- read.csv("../data/implied_costs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

# Set parameters
T <- nrow(econ_series)

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

discount_rate <- 0.05
discount_fac <- 1/(1+discount_rate)

##### Which method of generating the fe_eqm path is better? fe_eqm is used to generate the open access policy functions.

# Regression: This is what was estimated. It adjusts for measurement error.
# fe_eqm <- econ_coefs[1,2] + econ_coefs[2,2]*econ_series$r_s + econ_coefs[3,2]*econ_series$Ft_Ft # calculate the path of the OA eqm condition from the calibrated regression

# Observed: this is what was observed. it does not adjust for measurement error, but may fit the observed data better.
fe_eqm <- observed_time_series$risk[which(observed_time_series$year>2005)] # use the observed collision probability as the path of the OA eqm condition

#####

# The implied series has the "true" economic costs and returns implied by the measurement error equation. The procedure brings in future economic information from the calibration. Combine this with the observed economic series.
comb_econ_series <- merge(econ_series,implied_econ_series,by=c("year")) 

# Merge the economic and physical data. The final period gets truncated because the theory-adjustment formula uses up one year, so pad the final observation with the last value in the series. For p, this value is observed. For F, this value is just a copy of the newly-penultimate value. This carries some implications for the measurement errors in the final period.
physecon <- merge(observed_time_series,comb_econ_series,by=c("year","risk"))
p <- physecon$pi_t
p <- c(p,econ_series$tot_rev[length(econ_series$tot_rev)])
F <- physecon$F_hat
F <- c(F,F[length(F)])
# Normalize pi_1 := 1. This helps convergence to a fixed tolerance because the numbers are smaller to begin with. In the end, we undo the normalization to get things back into the right units.
norm_const <- p[1]
p <- p/norm_const
F <- F/norm_const

# Construct launch constraint. In each period, the number of possible launches are the largest number of launch attempts (failures+successes) observed until that period.
observed_launches <- observed_time_series$launch_successes[which(observed_time_series$year>=start_year)]+observed_time_series$launch_failures[which(observed_time_series$year>=start_year)]
launch_constraint <- cummax(observed_launches)
