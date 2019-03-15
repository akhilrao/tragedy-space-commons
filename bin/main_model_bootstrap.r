##### Script to generate main model bootstrap sensitivity analysis for "Tragedy of the Space Commons" paper.
###

source("final_script_bootstrap.r") # defaults to 30x30 grid, 50 bootstrap samples

main_sim <- read.csv("../data/main_simulation.csv")
bootstrap_sims <- read.csv("../data/bootstrapped_simulation.csv")

main_small <- data.frame(year=main_sim$year, 
						launches.oa=main_sim$launches.oa,
						satellites.oa=main_sim$satellites.oa,
						debris.oa=main_sim$debris.oa,
						collision_rate.oa=main_sim$collision_rate.oa,
						launches.opt=main_sim$launches.opt,
						satellites.opt=main_sim$satellites.opt,
						debris.opt=main_sim$debris.opt,
						collision_rate.opt=main_sim$collision_rate.opt,
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
						bootstrap_draw=bootstrap_sims$bootstrap_draw)						

main_small$opt_tax_path <- (main_small$collision_rate.oa/main_small$satellites.oa - main_small$collision_rate.opt/main_small$satellites.opt)*F_over_horizon*norm_const*1e+9/main_small$satellites.oa
bootstrap_small$opt_tax_path <- (bootstrap_small$collision_rate.oa/bootstrap_small$satellites.oa - bootstrap_small$collision_rate.opt/bootstrap_small$satellites.opt)*F_over_horizon*norm_const*1e+9/bootstrap_small$satellites.oa

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
						theme_minimal() + ggtitle(paste0("Residual bootstrap simulations, ",B," draws\nOpen access launch projections"))
m_bs_small_long_bootstrap_oasats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\nOA satellite projections"))
m_bs_small_long_bootstrap_oadebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\nOA debris projections"))
m_bs_small_long_bootstrap_oacoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.oa"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Open access collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\n\nOA collision rate projections"))

dev.new()
grid.arrange(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot, ncol=2)

png(width=600,height=600,filename="../images/bootstrapped_openaccess_simulation_plot.png")
grid.arrange(m_bs_small_long_bootstrap_oalaunch_plot, m_bs_small_long_bootstrap_oasats_plot, m_bs_small_long_bootstrap_oadebs_plot, m_bs_small_long_bootstrap_oacoll_plot, ncol=2)
dev.off()

m_bs_small_long_bootstrap_optlaunch_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="launches.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal launch rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + 
						ggtitle(paste0("Residual bootstrap simulations, ",B," draws\nOptimal launch projections"))
m_bs_small_long_bootstrap_optsats_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="satellites.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal satellite stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\nOptimal satellite projections"))
m_bs_small_long_bootstrap_optdebs_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="debris.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal debris stock") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\nOptimal debris projections"))
m_bs_small_long_bootstrap_optcoll_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="collision_rate.opt"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.9, alpha=0.45, color="gray") + 
						ylab("Optimal collision rate") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("\nOptimal collision rate projections"))

m_bs_small_long_bootstrap_opttax_plot <- ggplot(data=m_bs_small_long[which(m_bs_small_long$time=="opt_tax_path"),], aes(x=year)) + 
						geom_line(aes(y=bs.value, group=as.factor(bootstrap_draw)), size=0.91, alpha=0.45, color="gray") + 
						ylab("Optimal satellite tax (nominal USD)") +
						geom_line(aes(y=m.value, group=as.factor(bootstrap_draw)),size=1) +
						theme_minimal() + ggtitle(paste0("Optimal satellite tax projections\n",B, " residual bootstrap draws"))


dev.new()
m_bs_small_long_bootstrap_opttax_plot

dev.new()
grid.arrange(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot, ncol=2)

png(width=300,height=300,filename="../images/bootstrapped_optimal_sattax_projection.png")
m_bs_small_long_bootstrap_opttax_plot
dev.off()

png(width=600,height=600,filename="../images/bootstrapped_optimal_simulation_plot.png")
grid.arrange(m_bs_small_long_bootstrap_optlaunch_plot, m_bs_small_long_bootstrap_optsats_plot, m_bs_small_long_bootstrap_optdebs_plot, m_bs_small_long_bootstrap_optcoll_plot, ncol=2)
dev.off()
