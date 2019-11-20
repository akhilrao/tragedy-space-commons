##### Script to set calibrated economic and physical parameters based on observed data and theory

# Read in data and parameter estimates
## bootstrapped parameters
risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")

risk_cal_set <- read.csv("../data/bootstrapped_risk_eqn_coefs.csv")
deblom_cal_set <- read.csv("../data/bootstrapped_debris_lom_coefs.csv")
accepted_set <- which(risk_cal_set$SD>0)
risk_cal_accepted <- risk_cal_set[accepted_set,]
deblom_cal_accepted <- deblom_cal_set[accepted_set,]

## non-bootstrapped values
econ_coefs <- read.csv("../data/econ_series_coefs.csv")
implied_econ_series <- read.csv("../data/implied_costs.csv")
satlom_cal <- read.csv("../data/calibrated_satellite_lom_coefs.csv")
econ_series <- read.csv("../data/econ_series.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")
MS_proj_rev <- read.csv("../data/avg_econ_return.csv")
MS_proj_total <- read.csv("../data/avg_econ_total.csv")
commercial_sector_growth <- read.csv("../data/commercial_sector_growth.csv")

png(width=500,height=400,filename=paste0("../images/implied_v_observed_costs.png"))
ggplot(data=implied_econ_series, aes(x=year)) + 
	geom_line(aes(y=F_t),size=1) +
	geom_line(aes(y=F_hat),size=1, linetype="dashed") +
	ylab("Billion USD") + xlab("Year") + theme_bw() +
	ggtitle("Observed vs. implied aggregate launch costs") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) 
dev.off()

# construct the satellite sector growth projections. projection_end is set in final_script.r.
projection_start <- MS_proj_rev$Year[nrow(MS_proj_rev)]+1
revenue_mean_CAGR <- mean((MS_proj_rev[-1,2]/MS_proj_rev[-nrow(MS_proj_rev),2] - 1))
revenue_growth <- cumprod(rep(1+revenue_mean_CAGR,length=length(c(projection_start:projection_end))))
revenue_projection <- data.frame(Year=c(projection_start:projection_end),Revenues=(MS_proj_rev[nrow(MS_proj_rev),2]*revenue_growth))
total_mean_CAGR <- mean((MS_proj_total[-1,2]/MS_proj_total[-nrow(MS_proj_total),2] - 1))
total_growth <- cumprod(rep(1+total_mean_CAGR,length=length(c(projection_start:projection_end))))
total_projection <- data.frame(Year=c(projection_start:projection_end),Total=(MS_proj_total[nrow(MS_proj_total),2]*total_growth))
MS_proj_rev <- rbind(MS_proj_rev,revenue_projection)
MS_proj_total<- rbind(MS_proj_total,total_projection)

# Set parameters
selected_years <- which(observed_time_series$year>=start_year)
S0 <- observed_time_series$payloads_in_orbit[which(observed_time_series$year==start_year)]
D0 <- observed_time_series$debris[which(observed_time_series$year==start_year)]

# estimated parameterization
satlom_cal_names <- as.character(satlom_cal[,1])
satlom_cal <- data.frame(parameters=t(c(satlom_cal[,2])))
colnames(satlom_cal) <- satlom_cal_names
avg_sat_decay <- satlom_cal$payloads_in_orbit # corresponds to just over 30 years on orbit: on average 5 year mission time + 25 year post-mission disposal compliance. Value estimated from statistical model for satellite law of motion. 
# decay coefficient is currently represented as survival rate rather than decay rate.

risk_cal_names <- names(risk_cal)[2:3]
risk_cal <- data.frame(parameters=c(risk_cal[,2:3]))
colnames(risk_cal) <- risk_cal_names
aSS <- risk_cal$S2
aSD <- risk_cal$SD

# risk_cal_names <- as.character(risk_cal[,1])
# risk_cal <- data.frame(parameters=t(c(risk_cal[,2])))
# colnames(risk_cal) <- risk_cal_names
# aSS <- mean(risk_cal_accepted$S2)
# aSD <- mean(risk_cal_accepted$SD)

deblom_cal_names <- as.character(deblom_cal[,1])
deblom_cal <- data.frame(parameters=t(c(deblom_cal[,2])))
colnames(deblom_cal) <- deblom_cal_names
aDDbDD <- 0 # this should be zero, to be consistent with "no Kessler Syndrome allowed". the estimated parameter is small so it's unlikely that setting it to zero messes with the fit within the sample support. setting to zero is convenient because it prevents any numerical instabilities from causing a blowup in debris.
# bSS <- mean(deblom_cal_accepted$SSfrags)#deblom_cal$SSfrags
# bSD <- mean(deblom_cal_accepted$SDfrags)#deblom_cal$SDfrags
# d <- 1-mean(deblom_cal_accepted$debris) #1-deblom_cal$debris
# m <- mean(deblom_cal_accepted$launch_successes)#deblom_cal$launch_successes
# asat_coef <- mean(deblom_cal_accepted$num_destr_asat)#deblom_cal$num_destr_asat
bSS <- deblom_cal$SSfrags
bSD <- deblom_cal$SDfrags
d <- 1-deblom_cal$debris
m <- deblom_cal$launch_successes
asat_coef <- deblom_cal$num_destr_asat

discount_rate <- 0.05
discount_fac <- 1/(1+discount_rate)

# The implied series has the "true" economic costs and returns implied by the measurement error equation. The procedure brings in future economic information from the calibration. Combine this with the observed economic series.
comb_econ_series <- merge(econ_series,implied_econ_series,by=c("year")) 

# Merge the economic and physical data. The final period gets truncated because the theory-adjustment formula uses up one year.
#One option is to pad the final observation with the last value in the series. For p, this value is observed. For F, this value is just a copy of the newly-penultimate value. This carries some implications for the measurement errors in the final period.
physecon <- merge(observed_time_series,comb_econ_series,by=c("year"))
p <- physecon$pi_t
F <- physecon$F_hat

# Append data with average case revenue projection from Morgan Stanley.
econ_proj <- merge(MS_proj_rev,MS_proj_total,by=c("Year"))
econ_proj$Costs <- econ_proj$Total - econ_proj$Revenues
econ_proj[,c(2,3,4)] <- econ_proj[,c(2,3,4)]/1000 #the raw data are in units of millions of dollars. this brings them to units of billions of dollars.
p <- c(p,econ_proj$Revenues)
F <- c(F,econ_proj$Costs)

# Normalize pi_1 := 1. This helps convergence to a fixed tolerance because the numbers are smaller to begin with. In the end, we undo the normalization to get things back into the right units.
norm_const <- p[1]
p <- p/norm_const
F <- F/norm_const

##### fe_eqm is used to generate the open access policy functions.

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


#### Plot projections and other calibration data

# Launch constraint
lc_dfrm <- data.frame(year=seq(from=start_year,length.out=length(launch_constraint)),observed_constraint=c(obs_launch_constraint, rep(NA,length=(length(launch_constraint)-length(obs_launch_constraint)))),projected_constraint=launch_constraint)
lc_path_base <- ggplot(data=lc_dfrm,aes(x=year))
lc_path <- lc_path_base + geom_line(aes(y=projected_constraint),linetype="dashed",color="blue",size=0.8) +
		geom_line(aes(y=observed_constraint),size=1.2) +							
		geom_vline(xintercept=2015,size=1.2,linetype="dashed") +
		ylab("Cumulative maximum possible launches") + theme_minimal() +
		ggtitle("Observed and projected launch constraint")


png(width=300,height=300,filename="../images/linear_trend_launch_constraint.png")
lc_path
dev.off()

# Industry revenues and costs
econ_data_plot_series_orig <- econ_series[,1:3]
colnames(econ_data_plot_series_orig) <- c("Year","Revenues","Costs")
econ_data_plot_series_proj <- econ_proj[,c(1,2,4)]
econ_data_plot_series <- rbind(econ_data_plot_series_orig, econ_data_plot_series_proj)
econ_data_plot_series_long <- reshape(data=econ_data_plot_series, idvar="Year", varying=c("Revenues","Costs"), v.name=c("value"), times=c("Total industry revenues","Total industry costs"), direction="long", new.row.names=1:1000)
colnames(econ_data_plot_series_long) <- c("Year","Variable","Value")

revcost_plot <- ggplot(data=econ_data_plot_series_long, aes(x=Year, y=Value)) + 
				geom_line(aes(group=as.factor(Variable),linetype=as.factor(Variable)),size=1) +
				ggtitle("Aggregate commercial satellite costs and revenues") +
				labs(linetype="") +
				xlab("Year") +
				ylab("Billion USD")	+
				geom_vline(xintercept=econ_proj$Year[1],linetype="dashed",size=1,color="darkgray") +
				theme_minimal() +
				scale_linetype_discrete(labels=c("\n\nCosts:\nlaunch,\nmanufacturing,\nsupport\n","Revenues:\ntelecom,\nimaging"))	+
				theme(text=element_text(family="Helvetica",size=22),
					axis.text.x=element_text(family="Helvetica",size=22),
					axis.text.y=element_text(family="Helvetica",size=22),
					plot.title=element_text(family="Helvetica",size=22),
					legend.text=element_text(family="Helvetica",size=19))

png(width=800,height=400,filename="../images/industry_revcost_plot.png")
revcost_plot
dev.off()

# Growth in commercial relative to govt space
csg_long <- gather(commercial_sector_growth, variable, value, Commercial.Infrastructure.and.Support.Industries:Non.U.S..Government.Space.Budgets, factor_key=TRUE)

csg_base <- ggplot(data=csg_long[union(which(csg_long$Year==2005),which(csg_long$Year==2015)),],aes(as.factor(Year),value))

csg_plot <-	csg_base +
			geom_bar(aes(fill=variable), position="dodge", stat="identity" ) +
			labs(fill="") +
			ggtitle("Commercial space revenues and government space budgets") +
			ylab("Billion USD") +
			xlab("Year") +
			theme_minimal() +
			scale_fill_viridis(discrete=TRUE, labels=c("Comm'l Infrastructure\n& Support Industries\nRevenues\n", "Comm'l Space\nProducts & Services\nRevenues\n","US Govt\nSpace Budgets\n","Non-US Govt\nSpace Budgets")) +
			theme(text=element_text(family="Helvetica",size=22),
				axis.text.x=element_text(family="Helvetica",size=22),
				axis.text.y=element_text(family="Helvetica",size=22),
				plot.title=element_text(family="Helvetica",size=22),
				legend.text=element_text(family="Helvetica",size=19))

#####
# Main text and Extended Data figures
#####

# Main text figure 1
png(width=800,height=600,filename="../images/main_text_figure_1.png")
plot_grid(revcost_plot,csg_plot,labels=c("a","b"),align="v",axis="2",nrow=2,rel_widths=c(3/5,2/5),label_size=20)
dev.off()

# Extended Data Figure 6
lc_coef_table <- summary(lc_model)$coefficients[,1:2]

rownames(lc_coef_table) <- c("Intercept","Time trend")
lc_coef_table <- round(lc_coef_table,2)

lc_coef_table_plot <- ggtexttable(lc_coef_table, 
                        theme = ttheme("mOrange"))

png(width=450,height=400,filename=paste0("../images/extended_data_figure_6.png"))
plot_grid(lc_path, lc_coef_table_plot, align="h",labels=c("a","b"),axis="1",nrow=2,ncol=1,rel_heights=c(0.75,0.25))
dev.off()
