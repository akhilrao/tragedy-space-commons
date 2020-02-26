##### Script to generate main model bootstrap sensitivity analysis for "Tragedy of the Space Commons" paper.
###

if((file.exists(paste0("../data/bootstrapped_simulation.csv"))==FALSE)||(force_bootstrap_recalculation==1)) {
	source("final_script_bootstrap_calculations.r") # this actually calculates the bootstrap analysis. uses the same grid settings as the main model. won't run unless the file with bootstrap outputs doesn't exist, or unless forced to.
}

main_sim <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
bootstrap_sims <- read.csv("../data/bootstrapped_simulation.csv")

#main_sim <- main_sim[which(main_sim$start_time.opt==14|main_sim$start_time.opt==29),] # select only main_sim rows where optimal management begins in 2020 and 2035, to be comparable to the headline numbers of the paper.

main_small <- data.frame(year=main_sim$year, 
						launches.oa=main_sim$launches.oa,
						satellites.oa=main_sim$satellites.oa,
						debris.oa=main_sim$debris.oa,
						collision_rate.oa=main_sim$collision_rate.oa,
						launches.opt=main_sim$launches.opt,
						satellites.opt=main_sim$satellites.opt,
						debris.opt=main_sim$debris.opt,
						collision_rate.opt=main_sim$collision_rate.opt,
						costs=main_sim$costs.opt,
						NPV.oa=(main_sim$fleet_vfn_path.oa/main_sim$satellites.oa)*norm_const,
						NPV.opt=(main_sim$fleet_vfn_path.opt/main_sim$satellites.oa)*norm_const,
						start_time.opt=main_sim$start_time.opt,
						bootstrap_draw=0)

bootstrap_small <- data.frame(year=bootstrap_sims$year, 
						launches.oa=bootstrap_sims$launches.oa,
						satellites.oa=bootstrap_sims$satellites.oa,
						debris.oa=bootstrap_sims$debris.oa,
						collision_rate.oa=bootstrap_sims$collision_rate.oa,
						launches.opt=bootstrap_sims$launches.opt,
						satellites.opt=bootstrap_sims$satellites.opt,
						debris.opt=bootstrap_sims$debris.opt,
						collision_rate.opt=bootstrap_sims$collision_rate.opt,
						costs=bootstrap_sims$costs.opt,
						NPV.oa=(bootstrap_sims$fleet_vfn_path.oa/bootstrap_sims$satellites.oa)*norm_const,
						NPV.opt=(bootstrap_sims$fleet_vfn_path.opt/bootstrap_sims$satellites.oa)*norm_const,
						start_time.opt=bootstrap_sims$start_time.opt,
						bootstrap_draw=bootstrap_sims$bootstrap_draw)

# the variable naming convention is a little confusing at first: "NPV.welfare.gain" is the gains in a single year (2040) and simulation from having switched to optimal management in a selected year, e.g. 2020 (year 14) or 2035 (year 29). "npv_welfare_loss" is the difference in those gains from having switched in (say) 2035 instead of 2020, i.e. npv_welfare_loss in 2040 from switching in 2035 instead of 2020 = NPV.welfare.gain_2020 - NPV.welfare.gain_2035. since switching in 2020 instead of 2035 would yield additional gains to the tune of that npv_welfare_loss value, we consider the npv_welfare_loss value as the "gains in 2040 from switching in 2020 rather than 2035".
main_small$NPV.welfare.gain <- main_small$NPV.opt-main_small$NPV.oa
bootstrap_small$NPV.welfare.gain <- bootstrap_small$NPV.opt-bootstrap_small$NPV.oa

# calculate the price of anarchy ratios (same thing as welfare gain, just as a ratio -- need to undo the "per-satellite" terms factored in at lines 21-22 and 36-37)
main_small$NPV.PoA <- (main_small$NPV.opt/main_small$NPV.oa)#*(main_small$satellites.opt/main_small$satellites.oa)
bootstrap_small$NPV.PoA <- (bootstrap_small$NPV.opt/bootstrap_small$NPV.oa)#*(bootstrap_small$satellites.opt/bootstrap_small$satellites.oa)

# calculate optimal tax (OUF) paths
main_small$opt_tax_path <- (main_small$collision_rate.oa/main_small$satellites.oa - main_small$collision_rate.opt/main_small$satellites.opt)*main_small$costs*norm_const*1e+9/main_small$satellites.oa
bootstrap_small$opt_tax_path <- (bootstrap_small$collision_rate.oa/bootstrap_small$satellites.oa - bootstrap_small$collision_rate.opt/bootstrap_small$satellites.opt)*bootstrap_small$costs*norm_const*1e+9/bootstrap_small$satellites.oa

m_small_long <- reshape(main_small, idvar=c("year","bootstrap_draw","start_time.opt"), times=(colnames(main_small)[-c(1,which(colnames(main_small)=="bootstrap_draw"|colnames(main_small)=="start_time.opt"))]), direction="long", v.names="m.value", varying=(colnames(main_small)[-c(1,which(colnames(main_small)=="bootstrap_draw"|colnames(main_small)=="start_time.opt"))]))

bootstrap_small <- bootstrap_small[bootstrap_small$year>=2020,]
bs_small_long <- reshape(bootstrap_small, idvar=c("year","bootstrap_draw","start_time.opt"), times=(colnames(bootstrap_small)[-c(1,which(colnames(bootstrap_small)=="bootstrap_draw"|colnames(bootstrap_small)=="start_time.opt"))]), direction="long", v.names="bs.value", varying=(colnames(bootstrap_small)[-c(1,which(colnames(bootstrap_small)=="bootstrap_draw"|colnames(bootstrap_small)=="start_time.opt"))]))

rownames(m_small_long) <- NULL
rownames(bs_small_long) <- NULL

m_bs_small_long <- merge(m_small_long, bs_small_long, by=c("year","time","bootstrap_draw","start_time.opt"), suffixes=c(".main",".bs"), all=TRUE)
rownames(m_bs_small_long) <- NULL
head(m_bs_small_long)

# only keep runs beginning in 2006 to be comparable to ED fig 4
m_bs_small_long_06 <- m_bs_small_long[m_bs_small_long$start_time.opt==14,]
# only keep runs beginning in 2020 to be comparable to headline numbers
m_bs_small_long <- m_bs_small_long[m_bs_small_long$start_time.opt==14,]

# print time path of tax in main run
main_run_taxes <- intersect(which(m_bs_small_long$time=="opt_tax_path"),which(m_bs_small_long$bootstrap_draw==0))
tax_series <- m_bs_small_long$m.value[main_run_taxes]
message("Optimal OUF starts at ",round(tax_series[2],0)," USD in 2020 and escalates to ",round(tax_series[length(tax_series)],0)," USD in 2039, growing at an annualized rate of ", round(((tax_series[length(tax_series)]/tax_series[2])^(1/(length(tax_series)-1)) - 1)*100, 2),"% per year.")
	
# generate dataframe of value of switching in 2040 from main_small (headline number) and from bootstrap_small (distribution around headline number).

## for the main simulation results in main_small
m_coihist_df <- main_small[which(main_small$year==2040&main_small$start_time.opt>=14),c("year","start_time.opt","NPV.welfare.gain","NPV.PoA","bootstrap_draw")]
m_coihist_df <- ddply(m_coihist_df, .(year, bootstrap_draw), transform, npv_welfare_loss=(NPV.welfare.gain[which(start_time.opt==14)]-NPV.welfare.gain)/1000 )
m_coihist_df <- m_coihist_df[m_coihist_df$start_time.opt==14,]
## for the bootstrap draw results in bootstrap_small
bs_coihist_df <- bootstrap_small[which(bootstrap_small$year==2040&bootstrap_small$start_time.opt>=14),c("year","start_time.opt","NPV.welfare.gain","NPV.PoA","bootstrap_draw")]
bs_coihist_df <- ddply(bs_coihist_df, .(year,bootstrap_draw), transform, npv_welfare_loss=(NPV.welfare.gain[which(start_time.opt==14)] - NPV.welfare.gain)/1000 )
bs_coihist_df <- bs_coihist_df[bs_coihist_df$start_time.opt==14,]
## combine them for use later
coihist_df <- rbind(m_coihist_df,bs_coihist_df)

message("The middle 95% of model-predicted welfare gains in 2040 to optimal OUF implementation in 2020 fall in between ",round(quantile(coihist_df$NPV.PoA,probs=c(0.025,0.975)),2)[[1]], " and ",round(quantile(coihist_df$NPV.PoA,probs=c(0.025,0.975)),2)[[2]], " times the BAU NPV.")

##### Generate figures

m_bs_small_long_bootstrap_hist_plot <- ggplot(data=coihist_df) + 
						# geom_histogram(aes(x=npv_welfare_loss),fill="gray",bins=25) +
						# geom_vline(xintercept=coihist_df$npv_welfare_loss[coihist_df$bootstrap_draw==0]) +
						geom_histogram(aes(x=NPV.PoA),fill="gray",bins=125) +
						geom_vline(xintercept=coihist_df$NPV.PoA[coihist_df$bootstrap_draw==0]) +
						theme_bw() + ggtitle("Distribution of NPV gains from beginning optimal mgmt in 2020") + 
						#xlab("NPV gains in 2040 (nominal trillion USD)")	+
						xlab("Ratio of optimal NPV in 2040 to BAU NPV in 2040")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) )

middle_90 <- intersect(which(coihist_df$NPV.PoA<=quantile(coihist_df$NPV.PoA,probs=1)[[1]]),which(coihist_df$NPV.PoA>=quantile(coihist_df$NPV.PoA,probs=0)[[1]]))
coihist_df_mid <- coihist_df[middle_90,]
m_bs_small_long_bootstrap_hist_plot_mid90 <- ggplot(data=coihist_df_mid) + 
						geom_histogram(aes(x=NPV.PoA),fill="gray",bins=30) +
						geom_vline(xintercept=coihist_df$NPV.PoA[coihist_df$bootstrap_draw==0]) +
						theme_bw() + ggtitle("Distribution of NPV gains from beginning optimal mgmt in 2020") + 
						xlab("Ratio of optimal NPV in 2040 to BAU NPV in 2040")	+
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=20),
					axis.text.y=element_text(family="Helvetica",size=20),
					plot.title=element_text(family="Helvetica",size=20),
					legend.text=element_text(family="Helvetica",size=20) )

m_bs_small_long_bootstrap_oalaunch_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="launches.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access launch rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Open access launch projections")+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_oasats_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="satellites.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access satellite projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_oadebs_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="debris.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access debris projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_oacoll_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="collision_rate.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access collision rate projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

m_bs_small_long_bootstrap_optlaunch_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="launches.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal launch rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + 
						ggtitle("Optimal launch projections")+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_optsats_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="satellites.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal satellite projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_optdebs_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="debris.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal debris projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_optcoll_plot <- ggplot(data=m_bs_small_long_06[which(m_bs_small_long$time=="collision_rate.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal collision rate projections"))+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )

# tax and value function time paths 
# An OUF implemented in t+1 alters launch decisions in t. So to change behavior in 2020, the regulator announces an OUF in 2021. We plot the OUF from the year it begins changing behavior, i.e. we show the fee paid in 2021 as the 2020 value. "shifted_tax" accomplishes this.
m_bs_tax_shift <- m_bs_small_long[which(m_bs_small_long$time=="opt_tax_path"),]
m_bs_tax_shift <- ddply(m_bs_tax_shift, ~bootstrap_draw, transform, m_shifted_tax=c(m.value[-1],NA) )
m_bs_tax_shift <- ddply(m_bs_tax_shift, ~bootstrap_draw, transform, bs_shifted_tax=c(bs.value[-1],NA) )
m_bs_tax_shift <- m_bs_tax_shift[-which(m_bs_tax_shift$year==2041),]

m_bs_small_long_bootstrap_opttax_plot <- ggplot(data=m_bs_tax_shift, aes(x=year)) + 
						geom_line(aes(y=bs_shifted_tax, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						geom_line(aes(y=m_shifted_tax, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Optimal orbital-use fee projections") +
						scale_y_continuous(name="Optimal OUF (nominal $/sat)", labels = scales::comma)+
				theme(text=element_text(family="Helvetica",size=15),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) )
m_bs_small_long_bootstrap_oavalue_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="NPV.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Open access fleet NPV (nominal billion USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Open access fleet NPV projections")
m_bs_small_long_bootstrap_optvalue_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="NPV.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Optimal fleet NPV (nominal billion USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Optimal fleet NPV projections")

##### Main Text Figure

if(counterfactual=="none") {
	# MT figure 2
	png(width=1250,height=650,filename=paste0("../images/main_text_figure_2.png"))
	upper_row <- plot_grid(npv_welf_paths,coi_plot,labels=c("a","b"),align="h",axis="1",nrow=1,rel_widths=c(3/5,2/5),label_size=20)
	lower_row <- plot_grid(opt_tax_path,m_bs_small_long_bootstrap_hist_plot_mid90,labels=c("c","d"),align="h",axis="1",nrow=1,rel_widths=c(1/2,1/2),label_size=20)
	plot_grid(upper_row,lower_row,labels=c("",""),align="h",axis="1",nrow=2)
	dev.off()

	##### SI Figures

	# SI figure 4
	png(width=800,height=500,filename="../images/SI_fig_4.png")
	plot_grid(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot,align="h",labels=c("a","b","c","d"),axis="1",nrow=2,rel_widths=c(0.5,0.5),label_size=15)
	dev.off()

	# SI figure 5
	png(width=800,height=500,filename="../images/SI_fig_5.png")
	plot_grid(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot,align="h",labels=c("a","b","c","d"),axis="1",nrow=2,rel_widths=c(0.5,0.5),label_size=15)
	dev.off()

	# SI figure 6
	png(width=400,height=300,filename="../images/SI_fig_6.png")
	m_bs_small_long_bootstrap_opttax_plot
	dev.off()

}

