##### Script to generate main model figures for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Read in main model results, make additional outcome variables
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in main model results, make additional outcome variables
#############################################################################

OA_OPT <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))

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

shift_up <- function(x,optx) {
	output <- c(optx[1],x[1:(length(x)-1)])
	return(output)
}

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

#############################################################################
# 2.  Generate individual figures
#############################################################################

OA_fit_color <- "red"
OPT_fit_color <- "blue"
OA_OPT_fit_size <- 0.8
data_size <- 1

# Historical fit figures
OA_OPT_base_hist <- ggplot(data=OA_OPT[intersect(which(OA_OPT$year<=2015),which(OA_OPT$start_time.opt==0)),],aes(x=year))
OA_OPT_launch_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=launches.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launch_successes),size=data_size) +				
	ylab("Yearly launch rate") + theme_bw() +
	ggtitle("Simulated historical series\nLaunches")	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_sat_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=satellites.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),size=data_size) +
	ylab("LEO satellite stock") + theme_bw() +
	ggtitle(" \n Satellites")	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_deb_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=debris.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),size=data_size) +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_risk_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

# Projection figures: starting in 2006 only
OA_OPT_base_proj <- ggplot(data=OA_OPT[which(OA_OPT$start_time.opt==0),],aes(x=year))
OA_OPT_launch_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=launches.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launch_successes),size=data_size) +					
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ylab("Satellites launched") + xlab("Year") + theme_bw() +
	ggtitle("Launches")	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_sat_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=satellites.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + xlab("Year") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=15),
		axis.text.y=element_text(family="Helvetica",size=15),
		plot.title=element_text(family="Helvetica",size=15),
		legend.text=element_text(family="Helvetica",size=15) ) + 
	ylim(limits = c(0, max(OA_OPT$satellites.oa,OA_OPT$payloads_in_orbit)))
OA_OPT_deb_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=debris.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("Year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_risk_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dashed",color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("Year") + theme_bw()	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=15),
		axis.text.y=element_text(family="Helvetica",size=15),
		plot.title=element_text(family="Helvetica",size=15),
		legend.text=element_text(family="Helvetica",size=15) ) + 
	ylim(limits = c(0, max(OA_OPT$collision_rate.oa/OA_OPT$satellites.oa)))

# Projection figures: starting in multiple years
OA_OPT_base_proj_all <- ggplot(data=OA_OPT[which(OA_OPT$start_year==2006|OA_OPT$start_year==2015|OA_OPT$start_year==2020|OA_OPT$start_year==2030|OA_OPT$start_year==2035),], aes(x=year))
OA_OPT_launch_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=launches.opt, group=as.factor(start_year), color=as.factor(start_year)),linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa),color="dark gray",size=OA_OPT_fit_size) +
	geom_line(aes(y=launch_successes),size=data_size) +					
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_color_viridis(discrete=TRUE,guide=FALSE)	+
	ylab("Satellites launched") + theme_bw() +
	ggtitle("Launches")	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_sat_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=satellites.opt, group=as.factor(start_year), color=as.factor(start_year)),linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa),color="dark gray",size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	labs(color="Optimal mgmt\nstart year") +
	scale_color_viridis(discrete=TRUE,guide=FALSE)	+
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) ) + 
	ylim(limits = c(0, max(OA_OPT$satellites.oa)))
OA_OPT_deb_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=debris.opt, group=as.factor(start_year), color=as.factor(start_year)),linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa),color="dark gray",size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_color_viridis(discrete=TRUE,guide=FALSE)	+
	labs(color="Optimal mgmt\nstart year") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
OA_OPT_risk_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=collision_rate.opt/satellites.opt, group=as.factor(start_year), color=as.factor(start_year)),linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa),color="dark gray",size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_color_viridis(discrete=TRUE,guide=FALSE)	+
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) ) + 
	ylim(limits = c(0, max(OA_OPT$collision_rate.oa/OA_OPT$satellites.oa)))

# Prices of Anarchy and optimal tax path: starting in every year
risk_proj <- ggplot(data=OA_OPT,aes(x=year))
risk_comps <- risk_proj + 
	geom_line(aes(y=riskPoA,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Collision risk Price of Anarchy") + xlab("year") + theme_bw() +
	ggtitle("Improvement in satellite safety from optimal management") +
	scale_color_viridis(discrete=TRUE)	+
				theme(text=element_text(family="Helvetica",size=10),
					axis.text.x=element_text(family="Helvetica",size=10),
					axis.text.y=element_text(family="Helvetica",size=10),
					plot.title=element_text(family="Helvetica",size=10),
					legend.text=element_text(family="Helvetica",size=10) )
flow_welf_loss <- risk_proj + 
	geom_line(aes(y=flow_welfare_loss,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=0,linetype="dashed",color="blue") +
	ylab("Flow welfare gap (open access-optimal) \n(undiscounted $1b)") + xlab("year") + theme_bw() +
	ggtitle("Historical cost of open access and optimal tax path") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
				theme(text=element_text(family="Helvetica",size=10),
					axis.text.x=element_text(family="Helvetica",size=10),
					axis.text.y=element_text(family="Helvetica",size=10),
					plot.title=element_text(family="Helvetica",size=10),
					legend.text=element_text(family="Helvetica",size=10) )
npv_welf_loss <- risk_proj + 
	geom_line(aes(y=npv_welfare_loss,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	ylab("NPV welfare loss from open access ($1b)") + xlab("year") + theme_bw() +
	ggtitle("") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
				theme(text=element_text(family="Helvetica",size=10),
					axis.text.x=element_text(family="Helvetica",size=10),
					axis.text.y=element_text(family="Helvetica",size=10),
					plot.title=element_text(family="Helvetica",size=10),
					legend.text=element_text(family="Helvetica",size=10) )

# An OUF implemented in t+1 alters launch decisions in t. So to change behavior in 2020, the regulator announces an OUF in 2021. We plot the OUF from the year it begins changing behavior, i.e. we show the fee paid in 2021 as the 2020 value. "shifted_tax" in "OA_OPT_tax_shift" accomplishes this.
OA_OPT_tax_shift <- OA_OPT[which(OA_OPT$start_time.opt==14),]
#OA_OPT_tax_shift <- OA_OPT_tax_shift[-nrow(OA_OPT_tax_shift),]
OA_OPT_tax_shift$shifted_tax <- OA_OPT[which(OA_OPT$start_time.opt==14),]$opt_tax_path[-1]

opt_tax_path <- ggplot(data=OA_OPT_tax_shift,aes(x=year)) + 
	geom_line(aes(y=shifted_tax),size=data_size) + theme_bw() +
	labs(color="") +
	ylab("Optimal OUF ($/sat)") + xlab("Year") +
	ggtitle("Optimal OUF path") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) 

npv_poa_path <- risk_proj + 
	geom_line(aes(y=NPVPoA,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	guides(color=FALSE) +
	labs(color="") +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Fleet NPV\nimprovement factor") + xlab("year") + theme_bw() +
	ggtitle("Improvement in global satellite fleet NPV\nfrom optimal management") +
	scale_color_viridis(discrete=TRUE) +
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=15),
		axis.text.y=element_text(family="Helvetica",size=15),
		plot.title=element_text(family="Helvetica",size=15),
		legend.text=element_text(family="Helvetica",size=15) )
risk_poa_path <- risk_proj + 
	geom_line(aes(y=riskPoA,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Collision risk\nimprovement factor") + xlab("year") + theme_bw() +
	ggtitle("Improvement in satellite safety\nfrom optimal management") +
	#labs(color="Optimal mgmt\nstart year") +
	labs(color="") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

boa_base_dfrm <- OA_OPT[union(union(which(OA_OPT$year==2030),which(OA_OPT$year==2035)),which(OA_OPT$year==2040)),c("npv_welfare_gain","start_year","year")]		
coi_base_dfrm <- boa_base_dfrm[which(boa_base_dfrm$start_year>2010),]
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, npv_welfare_loss=(npv_welfare_gain[which(start_year==2020)]-npv_welfare_gain)/1000 )

coi_plot_cols <- c("2020" = paste0(viridis(7)[4]), "2025" = paste0(viridis(7)[5]), "2030" = paste0(viridis(7)[6]), "2035" = paste0(viridis(7)[7]))

coi_plot <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss)) +
			geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
			labs(fill="Optimal mgmt\nstart year") +
			ggtitle("Permanent orbit use value loss in 2040") +
			ylab("Forgone fleet NPV\n(nominal $1t)") +
			xlab("Year") +
			theme_bw() +
			scale_discrete_manual(values=coi_plot_cols, aesthetics = c("fill")) +
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) )

coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, npv_welfare_loss_pc=(npv_welfare_loss/npv_welfare_gain[which(start_year==2020)])*100 )
coi_plot_pc <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss_pc)) +
			geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
			labs(fill="Optimal mgmt\nstart year") +
			ggtitle("Permanent orbit\nuse value\nloss in 2040") +
			ylab("Forgone fleet NPV (percentage of 2020 optimal mgmt NPV)") +
			xlab("Year") +
			theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")) +
			theme_bw() +
			scale_discrete_manual(values=coi_plot_cols, aesthetics = c("fill")) +
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

boa_fin_dfrm <- OA_OPT[which(OA_OPT$year==2040),c("NPVPoA","start_year","year")]
row.names(boa_fin_dfrm) <- NULL
boa_plot_pc <- ggplot(data=boa_fin_dfrm[which(boa_fin_dfrm$start_year>=2020),],aes(start_year,NPVPoA)) +
				geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
				ylab("Value of the space industry with optimal management\n(multiple of BAU NPV in 2040)") +
				xlab("\nOptimal mgmt\nstart year") +
				theme_bw() +
				coord_flip() +
				guides(fill=FALSE) +
				theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
				scale_discrete_manual(values=coi_plot_cols, aesthetics = c("fill")) +
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

# divide by 1000 to get units of trillion USD
npvwelfpaths_long <- reshape2::melt(data=OA_OPT[,c("year","npv_oa_welfare","npv_opt_welfare","start_time.opt")], value.name="value", id=c("year","start_time.opt")) # recasts OA_OPT into a long dataframe for visualization in Main Text fig 2
npvwelfpaths_long <- npvwelfpaths_long[!duplicated(npvwelfpaths_long[,c("year","variable","value")]),] # removes open-access path rows which are duplicates but for different start_time.opt values
npvwelfpaths_long$start_time.opt[which(npvwelfpaths_long$variable=="npv_oa_welfare")] <- Inf
npvwelfpaths_long <- npvwelfpaths_long[which(npvwelfpaths_long$start_time.opt>=14),] # removes paths where optimal management begins before 2020
npv_welf_paths_plotbase <- ggplot(data=npvwelfpaths_long[npvwelfpaths_long$year>=2015,],aes(x=year)) # base for the plot
npvwelfpath_labs <- c(paste(c(opt_start_year[opt_start_year>=2020],"BAU (never)"),sep=",")) # label names
npvwelfpath_cols <- c("14" = paste0(viridis(7)[4]), "19" = paste0(viridis(7)[5]), "24" = paste0(viridis(7)[6]), "29" = paste0(viridis(7)[7]), "Inf" = "black") # label colors
npv_welf_paths <- npv_welf_paths_plotbase + 
	geom_line(aes(y=value/1000,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	theme(plot.margin = unit(c(0, 0, 0.5, 0), "cm")) +
	ylab("Fleet NPV (nominal $1t)") + xlab("Year") + theme_bw() +
	ggtitle("NPV under optimal management and BAU open access") +
	scale_discrete_manual(aesthetics=c("color"), labels=npvwelfpath_labs,values=npvwelfpath_cols)	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) )

risk_proj_nolt <- ggplot(data=OA_OPT[which(OA_OPT$start_year==2020),],aes(x=year))

opt_tax_path_solo <- risk_proj_nolt + 
	geom_line(aes(y=opt_tax_path),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_bw() +
	ggtitle("Optimal satellite tax path") +
	#scale_color_viridis(discrete=TRUE,labels=c(2020))	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) 

opt_dev_tax_path_solo <- risk_proj_nolt + 
	geom_line(aes(y=opt_dev_tax_path),size=data_size) +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_bw() +
	ggtitle("Optimal satellite tax path") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) + 
	ylim(limits = c(min(OA_OPT$opt_dev_tax_path[which(OA_OPT$start_time.opt==14)]), max(OA_OPT$opt_dev_tax_path[which(OA_OPT$start_time.opt==14)])))

opt_dev_tax_path_all <- risk_proj + 
	geom_line(aes(y=opt_dev_tax_path,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_bw() +
	labs(color="Optimal mgmt\nstart year") +
	ggtitle("Optimal satellite tax path") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) + 
	ylim(limits = c(min(OA_OPT$opt_dev_tax_path), max(OA_OPT$opt_dev_tax_path)))

risk_proj_20xx <- ggplot(data=OA_OPT[which(OA_OPT$start_year==2010|OA_OPT$start_year==2020|OA_OPT$start_year==2035),],aes(x=year))

opt_dev_tax_path_comp <- risk_proj_nolt + 
	geom_line(aes(y=opt_tax_path, group=as.factor(start_year), linetype=as.factor(start_year)),size=data_size, color="dark gray") +
	geom_line(aes(y=opt_dev_tax_path, group=as.factor(start_year), linetype=as.factor(start_year)),size=data_size) +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_bw() +
	ggtitle("Optimal satellite tax path comparison") +
	guides(linetype=FALSE)	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) )

opt_dev_tax_path_comp_all <- risk_proj_20xx + 
	geom_line(aes(y=opt_dev_tax_path,group=as.factor(start_year),linetype=as.factor(start_year)),size=data_size) +
	geom_line(aes(y=opt_tax_path,group=as.factor(start_year),linetype=as.factor(start_year)),size=data_size, color="dark gray") +
	guides(linetype=FALSE) +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_bw() +
#	labs(color="Optimal mgmt\nstart year") +
	ggtitle("Optimal satellite tax paths") +
#	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))	+
	theme(text=element_text(family="Helvetica",size=15),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) + 
	ylim(limits = c(0, max(OA_OPT$opt_tax_path)))

#############################################################################
# 3.  Main text and Extended Data figures
#############################################################################

# Main text figure 2 no longer made here, as it uses a bootstrap figure. produced in main_model_bootstrap.r instead.
png(width=850,height=450,filename=paste0("../images/main_text_figure_2.png"))
plot_grid(npv_welf_paths,coi_plot,labels=c("a","b"),align="h",axis="1",nrow=1,rel_widths=c(3.5/5,1.5/5),label_size=20)
dev.off()

# Extended data figure 4
png(width=800,height=600,filename=paste0("../images/extended_data_figure_4.png"))
plot_grid(OA_OPT_launch_proj,OA_OPT_sat_proj,OA_OPT_risk_proj,OA_OPT_deb_proj,align="h",axis="1",labels=c("a","b","c","d"),nrow=2,rel_widths=c(1/2,1/2),label_size=15)
dev.off()
