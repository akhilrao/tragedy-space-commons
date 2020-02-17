##### Script to generate figures for collision avoidance analysis
###
# Script flow:
# 1. Read in counterfactuals with collision avoidance scenario
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in alternate discount rate results, make additional outcome variables
#############################################################################

OA_OPT_main <- read.csv("../data/2006_7_starts_remfrac_0_remstart_2029_main_simulation.csv")


files <- list.files(path="../data/counterfactuals/military/", pattern="*.csv", full.names=FALSE, recursive=FALSE)
setwd("../data/counterfactuals/military/")

container <- list()
container[[1]] <- cbind(calc_tax_path(OA_OPT_main),mil_S=0)
for(i in 1:length(files)) {
	input <- read.csv(files[i+1])
	container[[i+1]] <- cbind(calc_tax_path(input),mil_S=mil_S_vary[i])
}

OA_OPT_mil_cf <- rbindlist(container)

#############################################################################
# 2.  Generate individual figures
#############################################################################

OA_OPT_mil_cf.small <- OA_OPT_mil_cf[,c("year","start_year","start_time.opt",
				"launches.oa","satellites.oa","debris.oa","collision_rate.oa",
				"launches.opt","satellites.opt","debris.opt","collision_rate.opt",
				"NPVPoA","opt_tax_path","mil_S")]

OA_OPT_mil_cf.small <- ddply(OA_OPT_mil_cf.small, .(start_year), transform, shifted_tax=c(opt_tax_path[-1],0))
OA_OPT_mil_cf.small <- OA_OPT_mil_cf.small[-which(OA_OPT_mil_cf.small$year==2041),]

OA_OPT_fit_size <- 0.8
data_size <- 1

OA_OPT_base_proj <- ggplot(data=OA_OPT_mil_cf.small[which(OA_OPT_mil_cf.small$start_time.opt==0),],aes(x=year))

(OA_OPT_launch_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=launches.opt,group=as.factor(mil_S),color=as.factor(mil_S)),size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa,group=as.factor(mil_S),color=as.factor(mil_S)),linetype="dashed",size=OA_OPT_fit_size) +	
	ylab("Satellites launched") + xlab("") + theme_bw() +
	ggtitle("Launches")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ) +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))

(OA_OPT_sat_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=satellites.opt,group=as.factor(mil_S),color=as.factor(mil_S)),size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa,group=as.factor(mil_S),color=as.factor(mil_S)),linetype="dashed",size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + xlab("") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))


(OA_OPT_deb_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=debris.opt,group=as.factor(mil_S),color=as.factor(mil_S)),size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa,group=as.factor(mil_S),color=as.factor(mil_S)),linetype="dashed",size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("Year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ) +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))

(OA_OPT_risk_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=collision_rate.opt/satellites.opt,group=as.factor(mil_S),color=as.factor(mil_S)), size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa/satellites.oa,group=as.factor(mil_S),color=as.factor(mil_S)),linetype="dashed",size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("Year") + theme_bw()	+
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) )  +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))

(opt_tax_path <- OA_OPT_base_proj + 
	geom_line(aes(y=shifted_tax,group=as.factor(mil_S),color=as.factor(mil_S)),size=data_size) + theme_bw() +
	scale_y_continuous(name="Optimal OUF (nominal $/sat)", labels = scales::comma) +
	xlab("Year") +
	ggtitle("Optimal OUF path") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) )  +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))

(NPVPoA_path <- OA_OPT_base_proj + 
	geom_line(aes(y=NPVPoA,group=as.factor(mil_S),color=as.factor(mil_S)),size=data_size) + theme_bw() +
	labs(color="") +
	scale_y_continuous(name="Ratio of optimal NPV to BAU NPV", labels = scales::comma) +
	xlab("Year") +
	ggtitle("NPV gains from optimal management") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) )  +
	labs(color="Number of\nmilitary satellites") +
	scale_color_viridis(discrete=TRUE,labels=c(paste(c(0,mil_S_vary),sep=","))))

#############################################################################
# 3.  Generate figure panel
#############################################################################

png(width=1000,height=700,filename=paste0("../../../images/military_counterfactual.png"))
	plot_grid(OA_OPT_launch_proj,OA_OPT_sat_proj,OA_OPT_risk_proj,OA_OPT_deb_proj,opt_tax_path,NPVPoA_path,align="v",axis="1",labels=c("a","b","c","d","e","f"),nrow=2,ncol=3,rel_widths=c(1/3,1/3,1/3),label_size=25)
dev.off()
