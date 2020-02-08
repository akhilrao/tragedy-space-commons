##### Script to generate figures for collision avoidance analysis
###
# Script flow:
# 1. Read in counterfactuals with collision avoidance scenario
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in alternate discount rate results, make additional outcome variables
#############################################################################

OA_OPT_estimated <- read.csv("../data/2006_7_starts_remfrac_0_remstart_2029_main_simulation.csv")
OA_OPT_avoid <- read.csv(paste0("../data/counterfactuals/collision_avoidance/",opt_start_year[1],"_cf_avoidance_aSS_",round(log(aSS),1),"_aSD_",round(log(aSD),1),"_simulation.csv"))

OA_OPT_estimated <- calc_tax_path(OA_OPT_estimated)
OA_OPT_avoid <- calc_tax_path(OA_OPT_avoid)

setwd("../data/counterfactuals/collision_avoidance/")

#############################################################################
# 2.  Generate individual figures
#############################################################################

selection <- c("year","start_year","start_time.opt",
				"launches.oa","satellites.oa","debris.oa","collision_rate.oa",
				"launches.opt","satellites.opt","debris.opt","collision_rate.opt",
				"NPVPoA","opt_tax_path")
OA_OPT_est.small <- OA_OPT_estimated[,selection]
OA_OPT_avo.small <- OA_OPT_avoid[,selection]

OA_OPT_est.small <- ddply(OA_OPT_est.small, .(start_year), transform, shifted_tax=c(opt_tax_path[-1],0))
OA_OPT_est.small <- OA_OPT_est.small[-which(OA_OPT_est.small$year==2041),]

OA_OPT_avo.small <- ddply(OA_OPT_avo.small, .(start_year), transform, shifted_tax=c(opt_tax_path[-1],0))
OA_OPT_avo.small <- OA_OPT_avo.small[-which(OA_OPT_avo.small$year==2041),]

OA_OPT <- merge(OA_OPT_est.small,OA_OPT_avo.small,by=c("year","start_year","start_time.opt"), suffix=c(".est",".avo") )

OA_fit_color <- "red"
OPT_fit_color <- "blue"
OA_OPT_fit_size <- 0.8
data_size <- 1

OA_OPT_base_proj <- ggplot(data=OA_OPT[which(OA_OPT$start_time.opt==14),],aes(x=year))

(OA_OPT_launch_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=launches.opt.est),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.est),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +	
	geom_line(aes(y=launches.opt.avo),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.avo),color=OA_fit_color,size=OA_OPT_fit_size) +	
	ylab("Satellites launched") + xlab("") + theme_bw() +
	ggtitle("Launches")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ) )

(OA_OPT_sat_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=satellites.opt.est),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.est),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.opt.avo),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.avo),color=OA_fit_color,size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + xlab("") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ))

(OA_OPT_deb_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=debris.opt.est),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.est),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.opt.avo),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.avo),color=OA_fit_color,size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("Year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ))

(OA_OPT_risk_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=collision_rate.opt.est/satellites.opt.est),linetype="dashed",color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.est/satellites.oa.est),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.opt.avo/satellites.opt.avo),color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.avo/satellites.oa.avo),color=OA_fit_color,size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("Year") + theme_bw()	+
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) )

(opt_tax_path <- OA_OPT_base_proj + 
	geom_line(aes(y=shifted_tax.est),linetype="dashed",size=data_size) + theme_bw() +
	geom_line(aes(y=shifted_tax.avo),size=data_size) + theme_bw() +
	labs(color="") +
	scale_y_continuous(name="Optimal OUF (nominal $/sat)", labels = scales::comma) +
	xlab("Year") +
	ggtitle("Optimal OUF path") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) )

(NPVPoA_path <- OA_OPT_base_proj + 
	geom_line(aes(y=NPVPoA.est),linetype="dashed",size=data_size) + theme_bw() +
	geom_line(aes(y=NPVPoA.avo),size=data_size) + theme_bw() +
	labs(color="") +
	scale_y_continuous(name="Ratio of optimal NPV to BAU NPV", labels = scales::comma) +
	xlab("Year") +
	ggtitle("NPV gains from optimal management") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) )

#############################################################################
# 3.  Generate figure panel
#############################################################################

png(width=1000,height=700,filename=paste0("../../../images/collision_avoidance_counterfactual.png"))
	plot_grid(OA_OPT_launch_proj,OA_OPT_sat_proj,OA_OPT_risk_proj,OA_OPT_deb_proj,opt_tax_path,NPVPoA_path,align="v",axis="1",labels=c("a","b","c","d","e","f"),nrow=2,ncol=3,rel_widths=c(1/3,1/3,1/3),label_size=25)
dev.off()
