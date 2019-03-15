##### Script to generate main model figures for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Read in main model results, make additional outcome variables
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in main model results, make additional outcome variables
#############################################################################

OA_OPT <- read.csv("../data/main_simulation.csv")

# Price of Anarchy in terms of collision risk. 1 represents no loss to anarchy, larger numbers show larger losses from anarchy.
OA_OPT$riskPoA <- (OA_OPT$collision_rate.oa/OA_OPT$collision_rate.opt)*(OA_OPT$satellites.opt/OA_OPT$satellites.oa)
# Price of Anarchy in terms of flow welfare. 1 : no present gains or losses to anarchy, >1 : present losses to anarchy, <1 : present gains to anarchy.
OA_OPT$flowWelfPoA <- OA_OPT$fleet_flowv.opt/OA_OPT$fleet_flowv.oa 
# Price of Anarchy in terms of NPV of welfare. 1 : no permanent gains or losses to anarchy, >1 : permanent losses to anarchy, <1 : permanent gains to anarchy.
OA_OPT$NPVPoA <- OA_OPT$fleet_vfn_path.opt/OA_OPT$fleet_vfn_path.oa 

# Since we're using aggregate data we need to divide by the number of satellites to get things into per-satellite units.
OA_OPT$flow_welfare_loss <- (OA_OPT$fleet_flowv.oa/OA_OPT$satellites.oa - OA_OPT$fleet_flowv.opt/OA_OPT$satellites.opt)*norm_const
OA_OPT$npv_welfare_loss <- (OA_OPT$fleet_vfn_path.oa/OA_OPT$satellites.oa - OA_OPT$fleet_vfn_path.opt/OA_OPT$satellites.opt)*norm_const
OA_OPT$npv_welfare_gain <- (OA_OPT$fleet_vfn_path.opt/OA_OPT$satellites.opt - OA_OPT$fleet_vfn_path.oa/OA_OPT$satellites.oa)*norm_const

F_over_horizon <- OA_OPT$costs.opt
OA_OPT$opt_tax_path <- (OA_OPT$collision_rate.oa/OA_OPT$satellites.oa - OA_OPT$collision_rate.opt/OA_OPT$satellites.opt)*F_over_horizon*norm_const*1e+9/OA_OPT$satellites.oa # 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

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
	ylab("Yearly launch rate") + theme_minimal() +
	ggtitle("Simulated historical series\nLaunches")
OA_OPT_sat_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=satellites.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),size=data_size) +
	ylab("LEO satellite stock") + theme_minimal() +
	ggtitle(" \n Satellites")
OA_OPT_deb_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=debris.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),size=data_size) +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_minimal()
OA_OPT_risk_hist <- OA_OPT_base_hist + 
	geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("year") + theme_minimal()

# Projection figures: starting in 2006 only
OA_OPT_base_proj <- ggplot(data=OA_OPT[which(OA_OPT$start_time.opt==0)),],aes(x=year))
OA_OPT_launch_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=launches.opt),linetype="dotted",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launch_successes),size=data_size) +					
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ylab("Satellites launched") + theme_minimal() +
	ggtitle("Simulated historical and projected series\nLaunches")
OA_OPT_sat_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=satellites.opt),linetype="dotted",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle(" \n Satellites") +
	ylab("LEO satellite stock") + theme_minimal()
OA_OPT_deb_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=debris.opt),linetype="dotted",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_minimal()
OA_OPT_risk_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=collision_rate.opt/satellites.opt),linetype="dotted",color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("year") + theme_minimal()


# # Projection figures: starting in 2006 only
# OA_OPT_base_proj <- ggplot(data=OA_OPT[which(OA_OPT$start_time.opt==0)),],aes(x=year))
# OA_OPT_launch_proj <- OA_OPT_base_proj + 
# 	geom_line(aes(y=launches.opt, group=as.factor(start_time.opt), color=start_time.opt),linetype="dotted",size=OA_OPT_fit_size) +
# 	geom_line(aes(y=launches.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=launch_successes),size=data_size) +					
# 	geom_vline(xintercept=2015,size=1,linetype="dashed") +
# 	ylab("Satellites launched") + theme_minimal() +
# 	ggtitle("Simulated historical and projected series\nLaunches")
# OA_OPT_sat_proj <- OA_OPT_base_proj + 
# 	geom_line(aes(y=satellites.opt, group=as.factor(start_time.opt)),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=payloads_in_orbit),size=data_size) +
# 	geom_vline(xintercept=2015,size=1,linetype="dashed") +
# 	ggtitle(" \n Satellites") +
# 	ylab("LEO satellite stock") + theme_minimal()
# OA_OPT_deb_proj <- OA_OPT_base_proj + 
# 	geom_line(aes(y=debris.opt, group=as.factor(start_time.opt)),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=debris.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=debris),size=data_size) +
# 	geom_vline(xintercept=2015,size=1,linetype="dashed") +
# 	ggtitle("Debris") +
# 	ylab("LEO debris stock") + xlab("year") + theme_minimal()
# OA_OPT_risk_proj <- OA_OPT_base_proj + 
# 	geom_line(aes(y=collision_rate.opt/satellites.opt, group=as.factor(start_time.opt)),linetype="dashed",color=OPT_fit_color, size=OA_OPT_fit_size) +
# 	geom_line(aes(y=collision_rate.oa/satellites.oa),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
# 	geom_line(aes(y=risk.x/payloads_in_orbit),size=data_size) +
# 	geom_vline(xintercept=2015,size=1,linetype="dashed") +
# 	ggtitle("Collision probability") +
# 	ylab("LEO collision risk") + xlab("year") + theme_minimal()

# Prices of Anarchy and optimal tax path: starting in 2006 and onwards
risk_proj <- ggplot(data=OA_OPT,aes(x=year))
risk_comps <- risk_proj + 
	geom_line(aes(y=riskPoA, group=as.factor(start_time.opt), color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	scale_colour_hue(labels=c(paste(opt_start_year,sep=","))) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Collision risk Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in satellite safety from optimal management")
flow_welf_loss <- risk_proj + 
	geom_line(aes(y=flow_welfare_loss, group=as.factor(start_time.opt), color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	scale_colour_hue(labels=c(paste(opt_start_year,sep=","))) +
	geom_hline(yintercept=0,linetype="dashed",color="blue") +
	ylab("Flow welfare gap (open access-optimal) \n(undiscounted $1b)") + xlab("year") + theme_minimal() +
	ggtitle("Historical cost of open access and optimal tax path")
npv_welf_loss <- risk_proj + 
	geom_line(aes(y=npv_welfare_loss, group=as.factor(start_time.opt), color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	scale_colour_hue(labels=c(paste(opt_start_year,sep=","))) +
	ylab("NPV welfare loss from open access ($1b)") + xlab("year") + theme_minimal() +
	ggtitle("")
opt_tax_path <- risk_proj + 
	geom_line(aes(y=opt_tax_path, group=as.factor(start_time.opt), color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	scale_colour_hue(labels=c(paste(opt_start_year,sep=","))) +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_minimal() +
	ggtitle("Optimal satellite tax path")
npv_poa_path <- risk_proj + 
	geom_line(aes(y=NPVPoA, group=as.factor(start_time.opt), color=as.factor(start_time.opt)),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	scale_colour_hue(labels=c(paste(opt_start_year,sep=","))) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management")

coi_base_dfrm <- OA_OPT[union(union(which(OA_OPT$year==2030),which(OA_OPT$year==2035)),which(OA_OPT$year==2040)),c("npv_welfare_gain","start_time.opt","year")]
colnames(coi_base_dfrm)[2] <- "start_time"
coi_base_dfrm$start_year <- rep(-1,length=nrow(coi_base_dfrm))
for(s in 1:length(unique(coi_base_dfrm$start_time))){
	coi_base_dfrm$start_year[coi_base_dfrm$start_time==sort(unique(coi_base_dfrm$start_time),decreasing=FALSE)[s]] <- unique(OA_OPT$year)[sort(unique(coi_base_dfrm$start_time),decreasing=FALSE)[s]+1]
}

coi_plot <- ggplot(data=coi_base_dfrm,aes(as.factor(year),npv_welfare_gain)) +
			geom_bar(aes(fill=as.factor(start_year), color=as.factor(start_year)), position="dodge", stat="identity" ) +
			scale_colour_hue(guide=FALSE) +
			labs(fill="Optimal mgmt\nstart year") +
			ggtitle("Costs of inaction:\nValue forgone from delaying start of optimal mgmt") +
			ylab("Lost social NPV") +
			xlab("Year") +
			theme_bw()			

#############################################################################
# 3.  Generate figure panels
#############################################################################

png(width=450,height=600,filename=paste0("../images/welfare_and_tax_optstart_",opt_start_year,"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
grid.arrange(opt_tax_path, risk_comps, npv_poa_path, ncol=1)
dev.off()


# historical fit path panel
png(width=500,height=500,filename=paste0("../images/simulated_historical_series_optstart_",opt_start_year,"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
grid.arrange(OA_OPT_launch_hist,OA_OPT_sat_hist,OA_OPT_risk_hist,OA_OPT_deb_hist,ncol=2)
dev.off()

# projected fit panel
png(width=500,height=500,filename=paste0("../images/simulated_projected_series_optstart_",opt_start_year,"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
grid.arrange(OA_OPT_launch_proj,OA_OPT_sat_proj,OA_OPT_risk_proj,OA_OPT_deb_proj,ncol=2)
dev.off()
