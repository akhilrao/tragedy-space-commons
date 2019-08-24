##### Script to generate main model figures for "Tragedy of the Space Commons" paper.
###
# Script flow:
# 1. Read in main model results, make additional outcome variables
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in main model results, make additional outcome variables
#############################################################################

R_frac <- D_fraction_to_remove

OA_OPT_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
OA_OPT_no_removal <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_0_remstart_",R_start_year,"_main_simulation.csv"))

#####

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
# 1e+9 scales to units of billion (nominal) dollars. "norm_const" is the normalization constant used during calibration to rescale the economic parameters (the rescaling makes the value function iteration better-behaved). We divide by the number of satellites to get the rate into a probability. The final division by the number of open access satellites converts the cost (F_over_horizon*norm_const*1e+9) from total dollars paid by industry into dollars per open-access satellite.

ss_rows <- which(OA_OPT$start_time.opt==-1)

OA_OPT$start_year <- rep(-1,length=nrow(OA_OPT))
OA_OPT_SS <- OA_OPT[ss_rows,c("year","npv_opt_welfare.rem","npv_opt_welfare.norem")]
colnames(OA_OPT_SS)[2:3] <- c("ss_npv_opt_welfare.rem","ss_npv_opt_welfare.norem")
OA_OPT <- OA_OPT[-ss_rows,]

# relabel start_time.opt values from integer labels to the appropriate years
for(s in 1:length(unique(OA_OPT$start_time.opt))){
	OA_OPT$start_year[OA_OPT$start_time.opt==sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]] <- unique(OA_OPT$year)[sort(unique(OA_OPT$start_time.opt),decreasing=FALSE)[s]+1]
}

OA_OPT <- merge(OA_OPT,OA_OPT_SS,by=c("year"))


# calculate costs of inaction (coi)
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

# npv_welfare_loss.rem is the cost of inaction (forgone gain from switching to opt mgmt in t instead of 2020, where t > 2020) with removal, .norem is without removal
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, npv_welfare_loss.norem=(npv_welfare_gain.norem[which(start_year==2020)]-npv_welfare_gain.norem), npv_welfare_loss.rem=(npv_welfare_gain.rem[which(start_year==2020)]-npv_welfare_gain.rem) )
# coi_effect_of_removal is the change in the cost of inaction due to introducing debris removal
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal=npv_welfare_loss.rem-npv_welfare_loss.norem)
# coi_effect_of_removal_pc is that change in cost of inaction expressed as a percentage of the gain from switching in 2020 with removal. "coi_effect_of_removal_pc=(coi_effect_of_removal/npv_welfare_gain.rem[which(start_year==2020)])*100" is the answer to, "supposing we had switched in 2020, with removal expected to come online in R_start_year, how much does the cost-of-inaction of switching in year t change when removal comes online?" so coi_effect_of_removal_pc = 7% says, "the introduction of debris removal in R_start_year increases the cost-of-inaction from switching in t by 7%, relative to having switched in 2020 with removal coming online in R_start_year"
coi_base_dfrm <- ddply(coi_base_dfrm, .(year), transform, coi_effect_of_removal_pc=(coi_effect_of_removal/npv_welfare_gain.rem[which(start_year==2020)])*100 )

#############################################################################
# 2.  Generate individual figures
#############################################################################

OA_fit_color <- "red"
OPT_fit_color <- "blue"
OA_OPT_fit_size <- 0.8
data_size <- 1

# Projection figures: starting in all years
OA_OPT_base_proj_all <- ggplot(data=OA_OPT[which(OA_OPT$start_year==2006),], aes(x=year))
OA_OPT_launch_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=launches.opt.rem), color="blue",linetype="dashed",size=(OA_OPT_fit_size)) +
	geom_line(aes(y=launches.oa.rem), color="red",linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.norem), color="red",size=OA_OPT_fit_size) +
	#geom_line(aes(y=launch_successes),color="dark gray",linetype="dotted",size=data_size) +
	geom_vline(xintercept=R_start_year,size=0.5,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ylab("Satellites launched") + theme_bw() +
	ggtitle("Launches")

OA_OPT_sat_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=satellites.opt.rem), color="blue",linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.rem),color="red",size=OA_OPT_fit_size,linetype="dashed") +
	geom_line(aes(y=satellites.oa.norem),color="red",size=OA_OPT_fit_size) +
	#geom_line(aes(y=payloads_in_orbit),color="dark gray",linetype="dotted",size=data_size) +
	geom_vline(xintercept=R_start_year,size=0.5,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + theme_bw()

OA_OPT_deb_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=debris.opt.rem, group=as.factor(start_year)), color="blue",linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.opt.norem, group=as.factor(start_year)), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.rem),color="red",size=OA_OPT_fit_size,linetype="dashed") +
	geom_line(aes(y=debris.oa.norem),color="red",size=OA_OPT_fit_size) +
	#geom_line(aes(y=debris),color="dark gray",linetype="dotted",size=data_size) +
	geom_vline(xintercept=R_start_year,size=0.5,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("year") + theme_bw()

OA_OPT_risk_proj_all <- OA_OPT_base_proj_all + 
	geom_line(aes(y=collision_rate.opt.rem/satellites.opt.rem), color="blue",linetype="dashed",size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.opt.norem/satellites.opt.norem), color="blue",size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.rem/satellites.oa.rem),color="red",size=OA_OPT_fit_size,linetype="dashed") +
	geom_line(aes(y=collision_rate.oa.norem/satellites.oa.norem),color="red",size=OA_OPT_fit_size) +
	#geom_line(aes(y=risk.x/payloads_in_orbit),color="dark gray",linetype="dotted",size=data_size) +
	geom_vline(xintercept=R_start_year,size=0.5,linetype="dashed") +
	scale_colour_hue(guide=FALSE) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Collision probability")+
	theme(axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20)) +
	ylab("LEO collision risk") + xlab("year") + theme_bw()

# cost-of-inaction (coi) figures
coi_plot.norem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss.norem)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Lost LEO\nuse value\nin 2040\n(no removal)") +
	ylab("Forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))
coi_plot.rem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),npv_welfare_loss.rem)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle("Lost LEO\nuse value\n in 2040\n(removal starting in 2030)") +
	ylab("Forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))

risk_proj <- ggplot(data=OA_OPT,aes(x=year))
npv_welf_paths <- risk_proj + 
	geom_line(aes(y=npv_oa_welfare.norem),size=data_size) +
	geom_line(aes(y=npv_oa_welfare.rem),size=data_size) +
	geom_line(aes(y=ss_npv_opt_welfare.norem),size=data_size,color=OPT_fit_color) +
	geom_line(aes(y=ss_npv_opt_welfare.rem),size=data_size*2,color=OPT_fit_color,alpha=0.6) +
	geom_line(aes(y=npv_opt_welfare.norem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size) +
	geom_line(aes(y=npv_opt_welfare.rem,group=as.factor(start_time.opt),color=as.factor(start_time.opt)),size=data_size*2,linetype="dashed") +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Social NPV (nominal $1b)") + xlab("Year") + theme_bw() +
	#ggtitle("Gains from shifting to optimal management\nrelative to open access BAU path (red line) and always-optimal path (black line)") +
	ggtitle("NPV gains of orbit recovery:\nshifting to optimal management from BAU open access\n(blue line: always-optimal, black line: BAU)") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(opt_start_year,sep=",")))
npv_oa_welf_paths <- risk_proj + 
	geom_line(aes(y=npv_oa_welfare.rem-npv_oa_welfare.norem),size=data_size) +
	labs(color="Optimal mgmt\nstart year") +
	ylab("Effect of debris removal on open access social NPV (nominal $1b)") + xlab("Year") + theme_bw() +
	ggtitle("NPV effect of debris removal on open access:")

coi_plot.deltarem <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),coi_effect_of_removal)) +
	geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
	labs(fill="Optimal mgmt\nstart year") +
	ggtitle(paste0("Change in lost LEO value for space industry in 2040\ndue to free removal beginning in ",R_start_year)) +
	ylab("Change in forgone social NPV (nominal $1b)") +
	xlab("Year") +
	theme_bw() +
	scale_fill_viridis(discrete=TRUE) +
	theme(legend.text=element_text(size=15))

coi_plot_cols <- c("2025"=viridis(3)[1],"2030"=viridis(3)[2],"2035"=viridis(3)[3])
coi_plot_pc <- ggplot(data=coi_base_dfrm[intersect(which(coi_base_dfrm$start_year>2020),which(coi_base_dfrm$year==2040)),],aes(as.factor(year),coi_effect_of_removal_pc)) +
			geom_bar(aes(fill=as.factor(start_year)), position="dodge", stat="identity" ) +
			labs(fill="Mgmt\nstart year") +
			ggtitle("Change in\nopen-access welfare loss \nin 2040") +
			ylab("Percentage of 2020 optimal mgmt NPV in 2040") +
			xlab("Year") +
			theme_bw() +
			scale_discrete_manual(values=coi_plot_cols, aesthetics = c("fill")) +
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

npv_poa_oa_remvnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.oaremvnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	labs(fill="Optimal mgmt\nstart year") +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_bw() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)
npv_poa_opt_remvnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.optremvnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_bw() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)

npv_poa_oaremvoptnorem_path <- risk_proj + 
	geom_line(aes(y=NPVPoA.oaremvoptnorem,group=as.factor(start_year),color=as.factor(start_year)),size=data_size) +
	guides(color=FALSE) +
	geom_hline(yintercept=1,linetype="dashed",color="blue") +
	ylab("Welfare NPV Price of Anarchy") + xlab("year") + theme_bw() +
	ggtitle("Improvement in global satellite fleet NPV from optimal management") +
	scale_color_viridis(discrete=TRUE)

#############################################################################
# 3.  Generate figure panels
#############################################################################

##### Nature submission figures

# Extended Data figure 5: removal projection fit panel, with benefits at the side

coi_total_summary <- read.csv(file="../data/pc_effect_of_removal_summary.csv")[,-1]
colnames(coi_total_summary) <- c("Largest percentage decrease\nin open-access welfare loss","Average percentage change\nin open-access welfare loss","Largest percentage increase\nin open-access welfare loss")
coi_total_summary <- round(coi_total_summary,2)#,"%")

removal_summary_table <- ggtexttable(coi_total_summary, rows = NULL, 
                        theme = ttheme("mOrange"))
png(width=800,height=500,filename=paste0("../images/extended_data_figure_5.png"))
ed5_left_column <- plot_grid(OA_OPT_launch_proj_all,OA_OPT_sat_proj_all,OA_OPT_risk_proj_all,OA_OPT_deb_proj_all,align="h",labels=c("a","b","c","d"),axis="1",nrow=2,rel_widths=c(0.5,0.5))
ed5_top_row <- plot_grid(ed5_left_column, coi_plot_pc, align="h",labels=c("","e"),axis="1",nrow=1,ncol=2,rel_widths=c(0.75,0.25))
plot_grid(ed5_top_row, removal_summary_table, align="h",labels=c("","f"),axis="1",nrow=2,ncol=1,rel_heights=c(0.85,0.15))
dev.off()
