B <- n_physical_bootstrap_draws # number of bootstrap resamples	

#####
# residual (semi-parametric) bootstrap for estimated nls parameters of physical collision risk model. grid searches same interval as original model estimate. 
#####
nls_risk_resids <- risk - nls_risk_xvars[,1]*(1 - exp(-nls_coefs[1,1]*nls_risk_xvars[,1] -nls_coefs[1,2]*nls_risk_xvars[,2]))

# generate bootstrap resamples of residuals 
resample_order <- replicate(B,sample(c(1:length(nls_risk_resids)),length(nls_risk_resids),replace = TRUE))	
nls_risk_resids_resamples <- matrix(nls_risk_resids[resample_order],ncol=B,byrow=FALSE)
nls_risk_model_output <- nls_risk_xvars[,1]*(1 - exp(-nls_coefs[1,1]*nls_risk_xvars[,1] -nls_coefs[1,2]*nls_risk_xvars[,2]))
nls_risk_bootstrap_dv <- nls_risk_model_output + nls_risk_resids_resamples

nls_risk_bootstrap_coefs <- matrix(-1,nrow=B,ncol=2)

# since the main-world estimates are the result of trying to find the best starting points for NLS, each set of bootstrap-world estimates should be too.
registerDoParallel(cores=ncores)
bootstrap_times <- rep(-1,length=B)
for(b in 1:B) {
	# search over a grid of initial points and select the best one
	message(paste0("Collision risk bootstrap loop ",b))
	
	bootstrap_clock_start <- proc.time()[3]
	
	## generate the grid
	init_risk_bootstrap_cond_element <- seq(from=1e-9,to=1e-6,length.out=15)
	init_risk_bootstrap_cond_grid <- expand.grid(init_risk_bootstrap_cond_element,init_risk_bootstrap_cond_element)
	nls_risk_bootstrap_values <- as.data.frame(matrix(1e+25,nrow=nrow(init_risk_bootstrap_cond_grid),ncol=4))
	colnames(nls_risk_bootstrap_values) <- c("obj_fn","SS_parm","SD_parm","conv")

	## run the search: estimate NLS at each starting point, choose the estimate which minimizes the objective
	nls_risk_bootstrap_values <- foreach(i=1:nrow(init_risk_bootstrap_cond_grid), .export=ls(), .inorder=TRUE, .combine=rbind) %dopar% {
		nls_risk_bootstrap_result <- invisible(spg(par=c(init_risk_bootstrap_cond_grid[i,1],init_risk_bootstrap_cond_grid[i,2]), fn=objective, rate=nls_risk_bootstrap_dv[,b], S=nls_risk_xvars$sats, D=nls_risk_xvars$debs, lower=0, upper=0.5,control=list(ftol=1e-15, trace=FALSE), quiet=TRUE)) 
		c(nls_risk_bootstrap_result$value, nls_risk_bootstrap_result$par[1], nls_risk_bootstrap_result$par[2], nls_risk_bootstrap_result$convergence)
	}
	colnames(nls_risk_bootstrap_values) <- c("obj_fn","SS_parm","SD_parm","conv")
	nls_risk_bootstrap_values <- data.frame(nls_risk_bootstrap_values)

	## select the best one
	best_init <- which.min(nls_risk_bootstrap_values$obj_fn)
	nls_risk_bootstrap_result <- nls_risk_bootstrap_values[best_init,2:3]
	nls_risk_bootstrap_coefs[b,] <- c(nls_risk_bootstrap_result[1,1],nls_risk_bootstrap_result[1,2])

	bootstrap_clock_end <- proc.time()[3]
	bootstrap_times[b] <- bootstrap_clock_end - bootstrap_clock_start

	if(b%%100==0) {
		avg_bootstrap_time <- round(mean(bootstrap_times[b:(b-99)]),3)
		message("Average time per bootstrap loop: ",avg_bootstrap_time," seconds")
	}

}
stopImplicitCluster()

nls_risk_bootstrap_coefs <- data.frame(nls_risk_bootstrap_coefs)
colnames(nls_risk_bootstrap_coefs) <- c("S2","SD")

# write out bootstrapped coefficients
write.csv(nls_risk_bootstrap_coefs,file="../data/bootstrapped_risk_eqn_coefs.csv")

#####
# residual (semi-parametric) bootstrap for estimated debris model parameters
#####
ridge_debris_resids <- D_next - m2xvars%*%m2coefs

# generate bootstrap resamples of residuals 
ridge_debris_resids_resamples <-  matrix(ridge_debris_resids[resample_order],ncol=B,byrow=FALSE)
ridge_debris_model_output <- m2xvars%*%m2coefs
ridge_debris_bootstrap_dv <- ridge_debris_model_output%*%matrix(1,ncol=B,nrow=1) + ridge_debris_resids_resamples

ridge_debris_bootstrap_coefs <- matrix(-1,nrow=B,ncol=ncol(m2xvars))

# create a cluster of workers estimating the debris model parameters on bootstrap-resampled data series
registerDoParallel(cores=ncores)
ridge_debris_bootstrap_coefs <- foreach(b=1:B, .export=ls(), .inorder=TRUE, .combine=rbind) %dopar% {
	bootstrap_dfrm <- dfrm #dfrm is the design matrix of physical variables from calibrate_physical_model.r
	Sfrac <- 1 - exp(-nls_risk_bootstrap_coefs[b,1]*bootstrap_dfrm$payloads_in_orbit)
	Dfrac <- 1 - exp(-nls_risk_bootstrap_coefs[b,2]*bootstrap_dfrm$debris)
	bootstrap_dfrm$SSfrags <- Sfrac*bootstrap_dfrm$payloads_in_orbit
	bootstrap_dfrm$SDfrags <- Dfrac*bootstrap_dfrm$payloads_in_orbit
	series <- cbind(bootstrap_dfrm[-(dim(bootstrap_dfrm)[1]),],D_next,S_next,risk)
	dfit_mat <- as.matrix(subset(series,select=c(payloads_in_orbit,payloads_decayed,debris,launch_successes,launch_failures,num_destr_asat,SSfrags,SDfrags) ) )
	dfit_mat_elnet <- as.matrix(subset(dfit_mat,select=c(debris,launch_successes,num_destr_asat,SSfrags,SDfrags) ) )
	m2_bootstrap <- glmnet(x=dfit_mat_elnet,y=D_next, alpha=0, lower.limits=c(0, 0, 0, 0, 0), upper.limits=c(1, Inf, Inf, Inf, Inf), lambda=cv.glmnet(x=dfit_mat_elnet,y=ridge_debris_bootstrap_dv[,b],alpha=0)$lambda.min,intercept=FALSE)
	coef(m2_bootstrap)[-1]
}
stopImplicitCluster()

ridge_debris_bootstrap_coefs <- data.frame(ridge_debris_bootstrap_coefs)
colnames(ridge_debris_bootstrap_coefs) <- colnames(m2xvars)
rownames(ridge_debris_bootstrap_coefs) <- NULL

# write out bootstrapped coefficients
write.csv(ridge_debris_bootstrap_coefs,file="../data/bootstrapped_debris_lom_coefs.csv")

#####
# Coefficient plots
#####

##### Collision risk equation
nls_risk_bootstrap_coefs <- read.csv("../data/bootstrapped_risk_eqn_coefs.csv")
positive_SDs <- which(nls_risk_bootstrap_coefs$SD>0)
effective_B <- length(positive_SDs)

# plot bootstrapped coefficient distribution
risk_hist_data <- nls_risk_bootstrap_coefs[positive_SDs,]
S2_hist <- ggplot(data=risk_hist_data, aes(x=S2)) + 
			geom_histogram(aes(y=..density..), bins=effective_B/3, colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
			xlab("Satellite-satellite collision parameter") +
			geom_vline(xintercept=nls_coefs[1,1], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(risk_hist_data$S2), linetype="dashed") +
			ggtitle(paste0("Distribution of bootstrapped satellite collision parameters: ",effective_B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
SD_hist <- ggplot(data=risk_hist_data, aes(x=SD)) + 
			geom_histogram(aes(y=..density..), bins=effective_B/3, colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
			xlab("Satellite-debris collision parameter") +
			geom_vline(xintercept=nls_coefs[1,2], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(risk_hist_data$SD), linetype="dashed") +
			ggtitle("\n") +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
grid.arrange(S2_hist,SD_hist,ncol=2)

png(width=525,height=550,filename="../images/collision_rate_parameter_bootstrap_plot.png")
grid.arrange(S2_hist,SD_hist,ncol=1)
dev.off()

nls_risk_bootstrap_dv_dfrm <- data.frame(nls_risk_bootstrap_dv)
colnames(nls_risk_bootstrap_dv_dfrm) <- paste0("bootstrapped_dv_",seq(from=1,to=B,by=1))

nls_risk_fitline <- nls_risk_model_output

nls_risk_fit <- data.frame(year=series$year, model_fit=nls_risk_fitline, nls_risk_bootstrap_dv)
nls_risk_fit_long <- reshape(nls_risk_fit, idvar="year", times=(colnames(nls_risk_fit)[-1]), direction="long", v.names="bootstrap_output", varying=(colnames(nls_risk_fit)[-1]))
nls_risk_fit_long$model_fit <- nls_risk_fit$model_fit
nls_risk_fit_long$truth <- risk
rownames(nls_risk_fit_long) <- NULL

##### Debris law of motion
ridge_debris_bootstrap_coefs <- read.csv("../data/bootstrapped_debris_lom_coefs.csv")

# condition on risk estimates both being positive
ridge_debris_bootstrap_coefs <- ridge_debris_bootstrap_coefs[positive_SDs,]

# draw histogram of coefficient distributions.
title_hist <- ggtitle("\n")
d_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=(1-debris))) + 
			geom_histogram(aes(y=..density..), bins=floor(effective_B/3), colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
			xlab("Debris decay") +
			geom_vline(xintercept=(1-m2coefs[1]), linetype="dashed", color="blue") +
			geom_vline(xintercept=(1-mean(ridge_debris_bootstrap_coefs$debris)), linetype="dashed") +
			ggtitle(paste0("Distribution of bootstrapped debris law of motion parameters: ", effective_B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
m_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=launch_successes)) + 
			geom_histogram(aes(y=..density..), bins=floor(effective_B/3), colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
			xlab("Launch debris") +
			geom_vline(xintercept=m2coefs[2], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$launch_successes), linetype="dashed") +
			ggtitle("\n") +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
gamma_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=num_destr_asat)) + 
			geom_histogram(aes(y=..density..), bins=floor(effective_B/3), colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
			xlab(paste0("Fragments from \nanti-satellite missile tests")) +
			geom_vline(xintercept=m2coefs[3], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$num_destr_asat), linetype="dashed") +
			ggtitle("\n") +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
beta_SS_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SSfrags)) + 
			geom_histogram(aes(y=..density..), bins=floor(effective_B/3), colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
			xlab(paste0("Fragments from \nsatellite-satellite collisions")) +
			geom_vline(xintercept=m2coefs[4], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SSfrags), linetype="dashed") +
			ggtitle("\n") +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )
beta_SD_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SDfrags)) + 
			geom_histogram(aes(y=..density..), bins=floor(effective_B/3), colour="gray", fill="white") +
			geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
			xlab(paste0("Fragments from \nsatellite-debris collisions")) +
			geom_vline(xintercept=m2coefs[5], linetype="dashed", color="blue") +
			geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SDfrags), linetype="dashed") +
			ggtitle("\n") +
			theme_bw()	+
			theme(text=element_text(size=15),
				axis.text.x=element_text(size=15),
				axis.text.y=element_text(size=15),
				plot.title=element_text(size=15),
				legend.text=element_text(size=15) )

# plot coefficient distributions
multiplot(title_hist, d_hist, m_hist, gamma_hist, beta_SS_hist, beta_SD_hist, 
			layout=matrix(c(2,1,3,4,5,6), 
			nrow=3, ncol=2, byrow=TRUE))
png(width=600,height=600,filename="../images/debris_bootstrap_parameter_plot.png")
multiplot(title_hist, d_hist, m_hist, gamma_hist, beta_SS_hist, beta_SD_hist, 
			layout=matrix(c(2,1,3,4,5,6), 
			nrow=3, ncol=2, byrow=TRUE))
dev.off()

ridge_debris_bootstrap_dv_dfrm <- data.frame(ridge_debris_bootstrap_dv)
colnames(ridge_debris_bootstrap_dv_dfrm) <- paste0("bootstrapped_dv_",seq(from=1,to=B,by=1))

ridge_debris_fitline <- ridge_debris_model_output

ridge_debris_fit <- data.frame(year=series$year, model_fit=ridge_debris_fitline, ridge_debris_bootstrap_dv)
ridge_debris_fit_long <- reshape(ridge_debris_fit, idvar="year", times=(colnames(ridge_debris_fit)[-1]), direction="long", v.names="bootstrap_output", varying=(colnames(ridge_debris_fit)[-1]))
ridge_debris_fit_long$model_fit <- ridge_debris_fit$model_fit
ridge_debris_fit_long$truth <- as.numeric(D_next)
rownames(ridge_debris_fit_long) <- NULL

#####
# Extended Data figures
#####

nls_risk_bootstrap_plot_base <- ggplot(data=nls_risk_fit_long, aes(x=year))
nls_risk_bootstrap_plot <- nls_risk_bootstrap_plot_base + geom_line(aes(y=bootstrap_output,group=as.factor(time)), size=0.5, alpha=0.2, color="light blue") + 
						geom_line(aes(y=model_fit), color="blue") +
						geom_line(aes(y=truth), size=1) +
						ylab("Expected number of collisions") +
						theme_minimal() + ggtitle("Collision rate calibration")

ridge_debris_bootstrap_plot_base <- ggplot(data=ridge_debris_fit_long, aes(x=year))
ridge_debris_bootstrap_plot <- ridge_debris_bootstrap_plot_base + geom_line(aes(y=bootstrap_output,group=as.factor(time)), size=0.5, alpha=0.2, color="light blue") + 
						ylab("Debris in LEO") +
						geom_line(aes(y=model_fit), color="blue") +
						geom_line(aes(y=truth), size=1) +
						theme_minimal() + ggtitle("Debris law of motion calibration")

# SI Figure 7

png(width=800,height=400,filename="../images/SI_figure_7.png")
plot_grid(nls_risk_bootstrap_plot,ridge_debris_bootstrap_plot,
		align="h",
		labels=c("a","b"),
		axis="1",
		nrow=1,
		rel_widths=c(1/2,1/2))
dev.off()
