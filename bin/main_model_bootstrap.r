##### Script to generate main model bootstrap sensitivity analysis for "Tragedy of the Space Commons" paper.
###

source("final_script_bootstrap_calculations.r") # this actually calculates the bootstrap analysis. uses the same grid settings as the main model.

main_sim <- read.csv(paste0("../data/",opt_start_year[1],"_",length(opt_start_year),"_starts_remfrac_",R_frac,"_remstart_",R_start_year,"_main_simulation.csv"))
bootstrap_sims <- read.csv("../data/bootstrapped_simulation.csv")

main_sim <- main_sim[which(main_sim$start_time.opt==0),] # select only main_sim rows where optimal management begins in 2006, to be comparable to the bootstrap draws.

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
						NPV.oa=main_sim$fleet_vfn_path.oa,
						NPV.opt=main_sim$fleet_vfn_path.opt,
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
						NPV.oa=bootstrap_sims$fleet_vfn_path.oa,
						NPV.opt=bootstrap_sims$fleet_vfn_path.opt,
						bootstrap_draw=bootstrap_sims$bootstrap_draw)						

main_small$opt_tax_path <- (main_small$collision_rate.oa/main_small$satellites.oa - main_small$collision_rate.opt/main_small$satellites.opt)*main_small$costs*norm_const*1e+9/main_small$satellites.oa
bootstrap_small$opt_tax_path <- (bootstrap_small$collision_rate.oa/bootstrap_small$satellites.oa - bootstrap_small$collision_rate.opt/bootstrap_small$satellites.opt)*bootstrap_small$costs*norm_const*1e+9/bootstrap_small$satellites.oa

m_small_long <- reshape(main_small, idvar=c("year","bootstrap_draw"), times=(colnames(main_small)[-c(1,ncol(main_small)-1)]), direction="long", v.names="m.value", varying=(colnames(main_small)[-c(1,ncol(main_small)-1) ]))
bs_small_long <- reshape(bootstrap_small, idvar=c("year","bootstrap_draw"), times=(colnames(bootstrap_small)[-c(1,ncol(bootstrap_small)-1)]), direction="long", v.names="bs.value", varying=(colnames(bootstrap_small)[-c(1,ncol(bootstrap_small)-1)]))
rownames(m_small_long) <- NULL
rownames(bs_small_long) <- NULL

m_bs_small_long <- merge(m_small_long, bs_small_long, by=c("year","time","bootstrap_draw"), suffixes=c(".main",".bs"), all=TRUE)
rownames(m_bs_small_long) <- NULL
head(m_bs_small_long)

m_bs_small_long_bootstrap_oalaunch_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="launches.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access launch rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Open access launch projections")
m_bs_small_long_bootstrap_oasats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access satellite projections"))
m_bs_small_long_bootstrap_oadebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access debris projections"))
m_bs_small_long_bootstrap_oacoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Open access collision rate projections"))

dev.new()
grid.arrange(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot, ncol=2)

png(width=600,height=600,filename="../images/bootstrapped_openaccess_simulation_plot.png")
grid.arrange(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot, ncol=2)
dev.off()

m_bs_small_long_bootstrap_optlaunch_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="launches.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal launch rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + 
						ggtitle("Optimal launch projections")
m_bs_small_long_bootstrap_optsats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal satellite projections"))
m_bs_small_long_bootstrap_optdebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal debris projections"))
m_bs_small_long_bootstrap_optcoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle(paste0("Optimal collision rate projections"))

m_bs_small_long_bootstrap_opttax_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="opt_tax_path"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Optimal satellite tax (nominal USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Optimal satellite tax projections")
m_bs_small_long_bootstrap_oavalue_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="NPV.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Open access fleet NPV (nominal USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Open access fleet NPV projections")
m_bs_small_long_bootstrap_optvalue_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="NPV.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Optimal fleet NPV (nominal USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_bw() + ggtitle("Optimal fleet NPV projections")

png(width=300,height=300,filename="../images/bootstrapped_optimal_sattax_projection.png")
m_bs_small_long_bootstrap_opttax_plot
dev.off()

png(width=600,height=600,filename="../images/bootstrapped_optimal_simulation_plot.png")
grid.arrange(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot, ncol=2)
dev.off()

##### Extended Data Figures

# ED figure 8
png(width=800,height=500,filename="../images/extended_data_figure_8.png")
plot_grid(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot,align="h",labels=c("a","b","c","d"),axis="1",nrow=2,rel_widths=c(0.5,0.5))
dev.off()

# ED figure 9
png(width=800,height=500,filename="../images/extended_data_figure_9.png")
plot_grid(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot,align="h",labels=c("a","b","c","d"),axis="1",nrow=2,rel_widths=c(0.5,0.5))
dev.off()

# ED figure 10
png(width=800,height=400,filename="../images/extended_data_figure_10.png")
plot_grid(m_bs_small_long_bootstrap_opttax_plot,m_bs_small_long_bootstrap_oavalue_plot,m_bs_small_long_bootstrap_optvalue_plot
,align="h",labels=c("a","b","c"),axis="1",nrow=1,rel_widths=c(1/3,1/3,1/3))
dev.off()

