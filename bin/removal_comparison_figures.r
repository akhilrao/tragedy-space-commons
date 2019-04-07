##### Script to generate main model figures for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Read in main model results, make additional outcome variables
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in main model results, make additional outcome variables
#############################################################################

OA_OPT_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
OA_OPT_no_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))

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


# Since we're using aggregate data we need to divide by the number of satellites to get things into per-satellite units.
OA_OPT$npv_oa_welfare.rem <- (OA_OPT$fleet_vfn_path.oa.rem/OA_OPT$satellites.oa.rem)*norm_const
OA_OPT$npv_opt_welfare.rem <- (OA_OPT$fleet_vfn_path.opt.rem/OA_OPT$satellites.opt.rem)*norm_const
OA_OPT$npv_welfare_loss.rem <- (OA_OPT$npv_oa_welfare.rem - OA_OPT$npv_opt_welfare.rem)
OA_OPT$npv_welfare_gain.rem <- (OA_OPT$npv_opt_welfare.rem - OA_OPT$npv_oa_welfare.rem)

OA_OPT$npv_oa_welfare.norem <- (OA_OPT$fleet_vfn_path.oa.norem/OA_OPT$satellites.oa.norem)*norm_const
OA_OPT$npv_opt_welfare.norem <- (OA_OPT$fleet_vfn_path.opt.norem/OA_OPT$satellites.opt.norem)*norm_const
OA_OPT$npv_welfare_loss.norem <- (OA_OPT$npv_oa_welfare.norem - OA_OPT$npv_opt_welfare.norem)
OA_OPT$npv_welfare_gain.norem <- (OA_OPT$npv_opt_welfare.norem - OA_OPT$npv_oa_welfare.norem)

OA_OPT$npv_oa_welfare_gain_from_rem <- (OA_OPT$fleet_vfn_path.oa.rem - OA_OPT$fleet_vfn_path.oa.norem)

F_over_horizon <- OA_OPT$costs.opt

OA_OPT$opt_tax_path.rem <- (OA_OPT$collision_rate.oa.rem/OA_OPT$satellites.oa.rem - OA_OPT$collision_rate.opt.rem/OA_OPT$satellites.opt.rem)*F_over_horizon*norm_const*1e+9 
OA_OPT$opt_tax_path.norem <- (OA_OPT$collision_rate.oa.norem/OA_OPT$satellites.oa.norem - OA_OPT$collision_rate.opt.norem/OA_OPT$satellites.opt.norem)*F_over_horizon*norm_const*1e+9 

# 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters for computational convenience. We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.


#############################################################################
# 2.  Generate individual figures
#############################################################################

OA_fit_color <- "red"
OPT_fit_color <- "blue"
OA_OPT_fit_size <- 0.8
data_size <- 1

ss_rows <- which(OA_OPT$start_time.opt==-1)

OA_OPT$start_year <- rep(-1,length=nrow(OA_OPT))
OA_OPT_SS <- OA_OPT[ss_rows,c("year","npv_opt_welfare.rem","npv_opt_welfare.norem")]
colnames(OA_OPT_SS)[2:3] <- c("ss_npv_opt_welfare.rem","ss_npv_opt_welfare.norem")
OA_OPT <- OA_OPT[-ss_rows,]

for(s in 1:length(unique(OA_OPT$start_time.opt))){
	OA_OPT$start_year[OA_OPT$start_time.opt==sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]] <- unique(OA_OPT$year)[sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]+1]
}

OA_OPT <- merge(OA_OPT,OA_OPT_SS,by=c("year"))


# Projection figures: starting in all years
OA_OPT_base_proj_all <- ggplot(data=OA_OPT[which(OA_OPT$start_year==2006),], aes(x=year))
OA_OPT_launch_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=launches.opt.rem), color="blue",linetype="dashed",size=(OA_OPT_fit_size*2)) +
	geom_line(aes(y=launches.oa.rem),linetype="dashed",size=OA_OPT_fit_size*2) +
	geom_line(aes(y=launches.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.norem),size=OA_OPT_fit_size) +
	geom_line(aes(y=launch_successes),color="dark gray",linetype="dotted",size=data_size*2) +					
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ylab("Satellites launched") + theme_minimal() +
	ggtitle("Simulated historical and projected series\nLaunches")

OA_OPT_sat_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=satellites.opt.rem), color="blue",linetype="dashed",size=OA_OPT_fit_size*2) +
	geom_line(aes(y=satellites.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.rem),size=OA_OPT_fit_size*2,linetype="dashed") +
	geom_line(aes(y=satellites.oa.norem),size=OA_OPT_fit_size) +
	geom_line(aes(y=payloads_in_orbit),color="dark gray",linetype="dotted",size=data_size*2) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle(" \n Satellites") +
	ylab("LEO satellite stock") + theme_minimal()

OA_OPT_deb_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=debris.opt.rem, group=as.factor(start_year)), color="blue",linetype="dashed",size=OA_OPT_fit_size*2) +
	geom_line(aes(y=debris.opt.norem, group=as.factor(start_year)), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.rem),size=OA_OPT_fit_size*2,linetype="dashed") +
	geom_line(aes(y=debris.oa.norem),size=OA_OPT_fit_size) +
	geom_line(aes(y=debris),color="dark gray",linetype="dotted",size=data_size*2) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_minimal()
OA_OPT_risk_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=collision_rate.opt.rem/satellites.opt.rem), color="blue",linetype="dashed",size=OA_OPT_fit_size*2) +
	geom_line(aes(y=collision_rate.opt.norem/satellites.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.rem/satellites.oa.rem),size=OA_OPT_fit_size*2,linetype="dashed") +
	geom_line(aes(y=collision_rate.oa.norem/satellites.oa.norem),size=OA_OPT_fit_size) +
	geom_line(aes(y=risk.x/payloads_in_orbit),color="dark gray",linetype="dotted",size=data_size*2) +
	geom_vline(xintercept=2015,size=1,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("year") + theme_minimal()

# Prices of Anarchy and optimal tax path: starting in every year
risk_proj <- ggplot(data=OA_OPT,aes(x=year))
risk_comps <- risk_proj + 
	geom_line(aes(y=riskPoA.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=riskPoA.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),linetype="dashed",size=(data_size*2)) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed") +
	ylab("Collision risk Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in satellite safety from optimal management") +
	scale_color_viridis(discrete=TRUE)
npv_welf_loss <- risk_proj + 
	geom_line(aes(y=npv_welfare_loss.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=npv_welfare_loss.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size*2,linetype="dashed") +
	labs(color="Optimal mgmt\nstart year") +
	ylab("NPV welfare loss from open access ($1b)") + xlab("year") + theme_minimal() +
	ggtitle("") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))
opt_tax_path <- risk_proj + 
	geom_line(aes(y=opt_tax_path.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=opt_tax_path.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size*2,linetype="dashed") +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Optimal satellite tax ($/sat)") + xlab("year") + theme_minimal() +
	ggtitle("Optimal satellite tax path") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))
npv_poa_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=NPVPoA.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size*2,linetype="dashed") +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)

#boa_base_dfrm <- OA_OPT[union(union(which(OA_OPT$year==2030),which(OA_OPT$year==2035)),which(OA_OPT$year==2040)),c("npv_welfare_gain.rem","npv_welfare_gain.norem","start_year","year")]
boa_base_dfrm <- OA_OPT[which(OA_OPT$year>=2030),c("npv_welfare_gain.rem","npv_welfare_gain.norem","start_year","year")]

boa_plot <- ggplot(data=boa_base_dfrm,aes(as.factor(year),npv_welfare_gain.rem)) +
			geom_bar(aes(fill=as.factor(start_year), color=as.factor(start_year)), position="dodge", stat="identity" ) +
			scale_colour_hue(guide=FALSE) +
			labs(fill="Optimal mgmt\nstart year") +
			ggtitle("Benefits of action:\nValue achieved by starting optimal mgmt early") +
			ylab("Gained social NPV") +
			xlab("Year") +
			theme_bw()			

coi_base_dfrm <- boa_base_dfrm[which(boa_base_dfrm$start_year>2010),]
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, npv_welfare_loss.norem=(npv_welfare_gain.norem[which(start_year==2015)]-npv_welfare_gain.norem), npv_welfare_loss.rem=(npv_welfare_gain.rem[which(start_year==2015)]-npv_welfare_gain.rem) )
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal=npv_welfare_loss.rem-npv_welfare_loss.norem)

coi_plot.norem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2015),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss.norem)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Lost LEO\nuse value\nin 2040\n(no removal)") +
	ylab("Forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))
coi_plot.rem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2015),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss.rem)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Lost LEO\nuse value\n in 2040\n(removal starting in 2030)") +
	ylab("Forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))
npv_welf_paths <- risk_proj + 
	geom_line(aes(y=npv_oa_welfare.norem),size=data_size) +
	geom_line(aes(y=npv_oa_welfare.rem),size=data_size) +
	geom_line(aes(y=ss_npv_opt_welfare.norem),size=data_size,color=OPT_fit_color) +
	geom_line(aes(y=ss_npv_opt_welfare.rem),size=data_size*2,color=OPT_fit_color,alpha=0.6) +
	geom_line(aes(y=npv_opt_welfare.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=npv_opt_welfare.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size*2,linetype="dashed") +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Social NPV (nominal $1b)") + xlab("Year") + theme_minimal() +
	#ggtitle("Gains from shifting to optimal management\nrelative to open access BAU path (red line) and always-optimal path (black line)") +
	ggtitle("NPV gains of orbit recovery:\nshifting to optimal management from BAU open access\n(blue line: always-optimal, black line: BAU)") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))
npv_oa_welf_paths <- risk_proj + 
	geom_line(aes(y=npv_oa_welfare.rem-npv_oa_welfare.norem),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Effect of debris removal on open access social NPV (nominal $1b)") + xlab("Year") + theme_minimal() +
	ggtitle("NPV effect of debris removal on open access:")

coi_plot.deltarem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2015),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),coi_effect_of_removal)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle(paste0("Change in lost LEO value for space industry in 2040\ndue to free removal beginning in ",R_start_year)) +
	ylab("Change in forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal_pc=(coi_effect_of_removal/npv_welfare_gain.rem[which(start_year==2015)])*100 )
coi_plot_pc <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2015),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),coi_effect_of_removal_pc)) +
			geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
			labs(fill="Optimal mgmt\nstart year") +
			ggtitle(paste0("Change in permanent orbit use value loss in 2040\ndue to removal beginning in ",R_start_year)) +
			ylab("Forgone fleet NPV (percentage of 2015 optimal mgmt NPV)") +
			xlab("Year") +
			theme_bw() +
			scale_discrete_manual(values=coi_plot_cols, aesthetics = c("fill")) +
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
min(coi_base_dfrm$coi_effect_of_removal_pc)
min(coi_base_dfrm$coi_effect_of_removal_pc[intersect(which(coi_base_dfrm$start_year>2015),which(coi_base_dfrm$year==2040))])

npv_poa_oa_remvnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.oaremvnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	labs(fill="Optimal mgmt\nstart year") +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)
npv_poa_opt_remvnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.optremvnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)

npv_poa_oaremvoptnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.oaremvoptnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_minimal() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)

#############################################################################
# 3.  Generate figure panels
#############################################################################

write.csv(coi_base_dfrm,file=paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_coi_base_dfrm.csv"))

png(width=500,height=500,filename=paste0("../images/remcomp_",length(opt_start_year),"_starts_change_in_welf_loss_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
coi_plot.deltarem
dev.off()

png(width=500,height=500,filename=paste0("../images/remcomp_",length(opt_start_year),"_starts_percent_change_in_welf_loss_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
coi_plot_pc
dev.off()

# policy benefits
png(width=600,height=750,filename=paste0("../images/remcomp_",length(opt_start_year),"_starts_welfare_and_tax_optstart_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
#grid.arrange(opt_tax_path, risk_comps, npv_poa_path, ncol=1)
plot_grid(opt_tax_path, risk_comps, npv_poa_path,align="h",axis="1",nrow=3,rel_heights=c(1/3,1/3,1/3))
dev.off()


# technology benefits
png(width=450,height=450,filename=paste0("../images/remcomp_",length(opt_start_year),"_starts_tech_benefits_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
#plot_grid(coi_plot.deltarem,npv_oa_welf_paths, align="h", axis="1", nrow=2,rel_heights=c(2,1))
plot_grid(npv_poa_oa_remvnorem_path,npv_poa_opt_remvnorem_path,align="v",axis="1",nrow=2,rel_heights=c(1,1))
dev.off()

# costs of inaction with removal
png(width=850,height=450,filename=paste0("../images/remcomp",length(opt_start_year),"_starts_costs_of_inaction_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
plot_grid(npv_welf_paths,coi_plot.norem,coi_plot.rem,npv_oa_welf_paths,NULL,NULL,align="h",axis="1",nrow=2,ncol=3,rel_widths=c(3/5,1/5,1/5),rel_heights=c(1/2,1/2))
dev.off()

# projected fit panel
png(width=600,height=600,filename=paste0("../images/remcomp",length(opt_start_year),"_starts_simulated_projected_series_optstart_",opt_start_year[1],"_remfrac_",R_frac,"_remstart_",R_start_year,".png"))
plot_grid(OA_OPT_launch_proj_all,OA_OPT_sat_proj_all,OA_OPT_risk_proj_all,OA_OPT_deb_proj_all,align="h",axis="1",nrow=2,rel_widths=c(1/2,1/2))
dev.off()
