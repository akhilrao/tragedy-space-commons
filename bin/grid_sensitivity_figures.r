##### Script to generate figures for grid resolution sensitivity analysis
###
# Script flow:
# 1. Read in outputs with different grid sizes
# 2. Modify and merge datasets
# 3. Calculate sensitivity values
# 4. Generate individual figures
# 5. Generate figure panels

#############################################################################
# 1.  Read in alternate grid size results, make additional outcome variables
#############################################################################

# Read in main model results
OA_OPT_coarse <- read.csv("../data/2006_7_starts_remfrac_0_remstart_2029_main_simulation_64pt.csv")
OA_OPT_fine <- read.csv("../data/2006_6_starts_remfrac_0_remstart_2029_main_simulation_900pt.csv")

# Read in bootstrap results
OA_OPT_bootstrap_coarse <- read.csv("../data/bootstrapped_simulation_64pt.csv")
OA_OPT_bootstrap_fine <- read.csv("../data/bootstrapped_simulation_900pt.csv")

# Attach optimal tax paths in case they're missing/NA

OA_OPT_coarse <- calc_tax_path(OA_OPT_coarse)
OA_OPT_fine <- calc_tax_path(OA_OPT_fine)
OA_OPT_bootstrap_coarse <- calc_tax_path(OA_OPT_bootstrap_coarse)
OA_OPT_bootstrap_fine <- calc_tax_path(OA_OPT_bootstrap_fine)

setwd("../data/counterfactuals/grid_size/")

#############################################################################
# 2.  Modify and merge datasets
#############################################################################

selection <- c("year","start_time.opt",
				"launches.oa","satellites.oa","debris.oa","collision_rate.oa",
				"launches.opt","satellites.opt","debris.opt","collision_rate.opt",
				"npv_opt_welfare", "npv_oa_welfare", "npv_welfare_gain",
				"NPVPoA","opt_tax_path")
OA_OPT_coarse.small <- OA_OPT_coarse[,selection]
OA_OPT_fine.small <- OA_OPT_fine[,selection]

OA_OPT_fine.small <- OA_OPT_fine.small %>% mutate(bootstrap_draw=0)
OA_OPT_coarse.small <- OA_OPT_coarse.small %>% mutate(bootstrap_draw=0)

OA_OPT_bs_coarse.small <- OA_OPT_bootstrap_coarse[,c(selection,"bootstrap_draw")]
OA_OPT_bs_fine.small <- OA_OPT_bootstrap_fine[,c(selection,"bootstrap_draw")]

# helper function to calculate shifted tax
shifted_tax_calc <- function(input) {
	output <- ddply(input, .(start_time.opt), transform, shifted_tax=c(opt_tax_path[-1],0))
	output <- output[-which(output$year==2041),]
}

OA_OPT_coarse.small <- shifted_tax_calc(OA_OPT_coarse.small)
OA_OPT_fine.small <- shifted_tax_calc(OA_OPT_fine.small)
OA_OPT_bs_coarse.small <- shifted_tax_calc(OA_OPT_bs_coarse.small)
OA_OPT_bs_fine.small <- shifted_tax_calc(OA_OPT_bs_fine.small)

OA_OPT_mbs.coarse <- merge(OA_OPT_coarse.small,OA_OPT_bs_coarse.small, by=c("year","bootstrap_draw","start_time.opt"), suffixes=c(".main",".bs"), all=TRUE)
OA_OPT_mbs.fine <- merge(OA_OPT_fine.small,OA_OPT_bs_fine.small, by=c("year","bootstrap_draw","start_time.opt"), suffixes=c(".main",".bs"), all=TRUE)

OA_OPT <- merge(OA_OPT_mbs.coarse,OA_OPT_mbs.fine,by=c("year","start_time.opt","bootstrap_draw"), suffix=c(".coarse",".fine") )

coihist_df_mbs.coarse <- rbind(OA_OPT_coarse.small,OA_OPT_bs_coarse.small)
coihist_df_mbs.fine <- rbind(OA_OPT_fine.small,OA_OPT_bs_fine.small)

coihist_df <- merge(coihist_df_mbs.coarse,coihist_df_mbs.fine,by=c("year","start_time.opt","bootstrap_draw"), suffix=c(".coarse",".fine") )

coihist_df <- coihist_df[which(coihist_df$year==2040&coihist_df$start_time.opt>=14) , c("year","start_time.opt",
	"NPVPoA.coarse", "NPVPoA.fine",
	"shifted_tax.coarse", "shifted_tax.fine",
	"bootstrap_draw")]

coihist_df <- coihist_df[coihist_df$start_time.opt==14,]

#############################################################################
# 3.  Calculate sensitivity values
#############################################################################

message("The main OUF estimate on the fine grid is ",round(coihist_df$shifted_tax.fine[which(coihist_df$bootstrap_draw==0)],2) ,", and on the coarse grid is ",round(coihist_df$shifted_tax.coarse[which(coihist_df$bootstrap_draw==0)],2))

message("The main optimal-to-BAU NPV ratio on the fine grid is ",round(coihist_df$NPVPoA.fine[which(coihist_df$bootstrap_draw==0)],2) ,", and on the coarse grid is ",round(coihist_df$NPVPoA.coarse[which(coihist_df$bootstrap_draw==0)],2))

message("The middle 95% of the optimal-to-BAU NPV ratio on the fine grid fall in between ",round(quantile(coihist_df$NPVPoA.fine,probs=c(0.025,0.975)),2)[[1]], " and ",round(quantile(coihist_df$NPVPoA.fine,probs=c(0.025,0.975)),2)[[2]], " times the BAU NPV.")

message("The middle 95% of the optimal-to-BAU NPV ratio on the coarse grid fall in between ",round(quantile(coihist_df$NPVPoA.coarse,probs=c(0.025,0.975)),2)[[1]], " and ",round(quantile(coihist_df$NPVPoA.coarse,probs=c(0.025,0.975)),2)[[2]], " times the BAU NPV.")


#############################################################################
# 4.  Generate individual figures
#############################################################################

OA_fit_color <- "red"
OPT_fit_color <- "blue"
OA_OPT_fit_size <- 0.8
data_size <- 1

OA_OPT_base_proj <- ggplot(data=OA_OPT[which(OA_OPT$bootstrap_draw==0&OA_OPT$start_time.opt==0),],aes(x=year))

(OA_OPT_launch_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=launches.opt.main.coarse),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.main.coarse),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +	
	geom_line(aes(y=launches.opt.main.fine),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=launches.oa.main.fine),color=OA_fit_color,size=OA_OPT_fit_size) +	
	ylab("Satellites launched") + xlab("") + theme_bw() +
	ggtitle("Launches")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ) )

(OA_OPT_sat_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=satellites.opt.main.coarse),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.main.coarse),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.opt.main.fine),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=satellites.oa.main.fine),color=OA_fit_color,size=OA_OPT_fit_size) +
	ggtitle("Satellites") +
	ylab("LEO satellite stock") + xlab("") + theme_bw() +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ))

(OA_OPT_deb_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=debris.opt.main.coarse),linetype="dashed",color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.main.coarse),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.opt.main.fine),color=OPT_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=debris.oa.main.fine),color=OA_fit_color,size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Debris") +
	ylab("LEO debris stock") + xlab("Year") + theme_bw()	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ))

(OA_OPT_risk_proj <- OA_OPT_base_proj + 
	geom_line(aes(y=collision_rate.opt.main.coarse/satellites.opt.main.coarse),linetype="dashed",color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.main.coarse/satellites.oa.main.coarse),linetype="dashed",color=OA_fit_color,size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.opt.main.fine/satellites.opt.main.fine),color=OPT_fit_color, size=OA_OPT_fit_size) +
	geom_line(aes(y=collision_rate.oa.main.fine/satellites.oa.main.fine),color=OA_fit_color,size=OA_OPT_fit_size) +
	#geom_vline(xintercept=2015,size=1,linetype="dashed") +
	ggtitle("Collision probability") +
	ylab("LEO collision risk") + xlab("Year") + theme_bw()	+
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) )

(opt_tax_path <- OA_OPT_base_proj + 
	geom_line(aes(y=shifted_tax.main.coarse),linetype="dashed",size=data_size) + theme_bw() +
	geom_line(aes(y=shifted_tax.main.fine),size=data_size) + theme_bw() +
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

(NPV_paths <- OA_OPT_base_proj + 
	geom_line(aes(y=npv_oa_welfare.main.coarse),linetype="dashed",color=OA_fit_color,size=data_size) + theme_bw() +
	geom_line(aes(y=npv_opt_welfare.main.coarse),linetype="dashed",color=OPT_fit_color,size=data_size) + theme_bw() +
	geom_line(aes(y=npv_oa_welfare.main.fine),color=OA_fit_color,size=data_size) + theme_bw() +
	geom_line(aes(y=npv_opt_welfare.main.fine),color=OPT_fit_color,size=data_size) + theme_bw() +
	labs(color="") +
	scale_y_continuous(name="Ratio of optimal NPV to BAU NPV", labels = scales::comma) +
	xlab("Year") +
	ggtitle("NPV paths on coarse and fine grids") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) )

(m_bs_small_long_bootstrap_hist_plot <- ggplot(data=coihist_df) +
						geom_histogram(aes(x=NPVPoA.fine),fill="darkgray",bins=50) +
						geom_vline(xintercept=coihist_df$NPVPoA.fine[coihist_df$bootstrap_draw==0]) +
						geom_histogram(aes(x=NPVPoA.coarse),fill="lightgray",bins=50) +
						geom_vline(xintercept=coihist_df$NPVPoA.coarse[coihist_df$bootstrap_draw==0], linetype="dashed", size=1) +
						theme_bw() + ggtitle("Distribution of NPV gains") + 
						#xlab("NPV gains in 2040 (nominal trillion USD)")	+
						xlab("Ratio of optimal NPV in 2040 to BAU NPV in 2040")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) ) )



#############################################################################
# 3.  Generate figure panel
#############################################################################

png(width=1400,height=800,filename=paste0("../../../images/SI_fig_11.png"))
	plot_grid(OA_OPT_launch_proj,OA_OPT_sat_proj,OA_OPT_risk_proj,OA_OPT_deb_proj,opt_tax_path,m_bs_small_long_bootstrap_hist_plot,align="v",axis="1",labels=c("a","b","c","d","e","f"),nrow=2,ncol=3,rel_widths=c(1/3,1/3,1/3),label_size=25)
dev.off()
