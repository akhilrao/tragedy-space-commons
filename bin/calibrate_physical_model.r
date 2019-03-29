#####################################################################
# Script to take a series of orbital stock levels and calibrate a physical model of debris evolution and collision risk
#####################################################################
#rm(list=ls())

system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

library(ggplot2)
library(glmnet)
library(grid)
library(gridExtra)
library(reshape2)
library(BB)
library(doParallel)

setwd("../data/")
set.seed(501) # for bootstrap samples

dfrm <- read.csv("ST_stock_series.csv")
# limit to years after 1990
#dfrm <- dfrm[which(dfrm$year>=1990),]
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]

##### define plotting functions
fitplot <- function(xvars,coefs,year,truth,title,ylabel) {
	fitline <- xvars%*%coefs
	fit <- data.frame(year=year,fit=fitline,truth=truth,error=(fitline-truth))

	plot_base <- ggplot(data=fit, aes(x=year))
	plot_fitplot <- plot_base + geom_point(aes(y=truth),size=1.1) +
							geom_line(aes(y=fit),size=0.9,linetype="dashed", color="blue") +
							theme_minimal() + ggtitle(paste(title)) +
							ylab(paste0(ylabel))	+
				theme(text=element_text(size=15),
					axis.text.x=element_text(size=15),
					axis.text.y=element_text(size=15),
					plot.title=element_text(size=15) )
	plot_error <- plot_base + geom_line(aes(y=error),size=0.9) +
						geom_hline(yintercept=0,linetype="dashed") +
						theme_minimal()	+
				theme(text=element_text(size=15),
					axis.text.x=element_text(size=15),
					axis.text.y=element_text(size=15),
					plot.title=element_text(size=15) )

	grid.arrange(plot_fitplot,plot_error,nrow=2)
}

nls_fitplot <- function(xvars,coefs,year,truth,title,ylabel) {
	fitline <- xvars[,1]*(1 - exp(-coefs[1]*xvars[,1] -coefs[2]*xvars[,2]))
	fit <- data.frame(year=year,fit=fitline,truth=truth,error=(fitline-truth))

	plot_base <- ggplot(data=fit, aes(x=year))
	plot_fitplot <- plot_base + geom_point(aes(y=truth),size=1.1) +
							geom_line(aes(y=fit),size=0.9,linetype="dashed", color="blue") +
							theme_minimal() + ggtitle(paste(title)) +
							ylab(paste0(ylabel))	+
				theme(text=element_text(size=15),
					axis.text.x=element_text(size=15),
					axis.text.y=element_text(size=15),
					plot.title=element_text(size=15) )
	plot_error <- plot_base + geom_line(aes(y=error),size=0.9) +
						geom_hline(yintercept=0,linetype="dashed") +
						theme_minimal()	+
				theme(text=element_text(size=15),
					axis.text.x=element_text(size=15),
					axis.text.y=element_text(size=15),
					plot.title=element_text(size=15) )

	grid.arrange(plot_fitplot,plot_error,nrow=2)
}

# Multiple plot function: taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##### Create additional X variables for regressions
dfrm$D2 <- dfrm$debris^2
dfrm$SD <- dfrm$payloads_in_orbit*dfrm$debris
dfrm$S2 <- dfrm$payloads_in_orbit^2

##### Create outcome variables for regressions
D_next <- c(dfrm$debris[-1])
S_next <- c(dfrm$payloads_in_orbit[-1])

##### Bind created variables into dataframe/matrix for regressions
series <- cbind(dfrm[-(dim(dfrm)[1]),],D_next,S_next)
ellfit_mat <- as.matrix(subset(series,select=c(S2,SD)))
ellfit_dfrm <- subset(series,select=c(S2,SD))
sfit_mat <- as.matrix(subset(series,select=c(payloads_in_orbit,payloads_decayed,launch_successes) ))

##### Run regressions to calibrate things

### satellite equation calibration. this is a sanity check: all coefficients should be +1 or -1 to indicate accounting relationship.
sfit_dfrm <- subset(series,select=c(payloads_in_orbit,payloads_decayed,launch_successes) )
sat_ols <- lm(S_next ~ -1 + payloads_in_orbit + payloads_decayed + launch_successes, data=sfit_dfrm)
summary(sat_ols)
sat_parms <- coef(sat_ols)

# plot fit
sat_xvars <- as.matrix(cbind(subset(series,select=colnames(sfit_mat))))
sat_coefs <- sat_parms
fitplot(sat_xvars,sat_coefs,series$year,S_next,"Satellite accounting model","satellites")
dev.off()

### satellite equation calibration. this is the statistical equation: calculating decay rate implied by including collisions. Produces a positive, large, and insignificant coefficient on risk, in addition to about 0.5 and significant coefficient on launch_successes. Implied satellite time on-orbit is close to 30 years, though, which seems reasonable.
# second approach (shown here): subtract off risk and launches, then calculate decay coefficient.
sfit_dfrm <- subset(series,select=c(payloads_in_orbit,risk,launch_successes) )

sat_ols <- lm(S_next ~ -1 + payloads_in_orbit + risk + launch_successes, data=sfit_dfrm)
summary(sat_ols)
sat_parms <- coef(sat_ols)

# plot fit
sat_xvars <- as.matrix(cbind(subset(sfit_dfrm,select=colnames(sfit_dfrm))))
sat_coefs <- c(sat_parms[1],-1*sat_parms[1],1)
fitplot(sat_xvars,sat_coefs,series$year,S_next,"Satellite statistical model (adjusted coefficients)","satellites")
dev.off()

fitplot(sat_xvars,sat_parms,series$year,S_next,"Satellite statistical model (unadjusted coefficients)","satellites")
dev.off()

write.csv(sat_parms[1],file="calibrated_satellite_lom_coefs.csv")

### satellite equation calibration. this is the physically-calibrated equation: decay rate set to about 30 years on orbit (from statistical approach above), losses and launches are 1:1.
# sat_coefs <- c(1-(1/sat_parms[1]), -1, 1)

### risk equation calibration. estimate coefficients by NLS.

# create outcome variable
risk <- series$risk

# NLS objective function for negative exponential collision rate
objective <- function(parms,rate,S,D) {
	aSS <- parms[1]
	aSD <- parms[2]
	obj <- sum((rate - S*(1 - exp(-(aSS*S)-(aSD*D))))^2)
	return(obj)
}

nls_risk_xvars <- subset(series,select=c(payloads_in_orbit,debris))
colnames(nls_risk_xvars) <- c("sats","debs")
print(find_best_nls_parms)
## Run through grid of initial guesses to get the best one if find_best_nls_parms = 1.
if(find_best_nls_parms == 1) {
	registerDoParallel(cores=ncores)
	init_cond_element <- seq(from=1e-9,to=1e-6,length.out=400)
	init_cond_grid <- expand.grid(init_cond_element,init_cond_element)
	nls_values <- as.data.frame(matrix(1e+25,nrow=nrow(init_cond_grid),ncol=4))
	colnames(nls_values) <- c("obj_fn","SS_parm","SD_parm","conv")

	nls_values <- foreach(i=1:nrow(init_cond_grid), .export=ls(), .inorder=TRUE, .combine=rbind) %dopar% {
		nls_result <- invisible(spg(par=c(init_cond_grid[i,1],init_cond_grid[i,2]), fn=objective, rate=risk, S=nls_risk_xvars$sats, D=nls_risk_xvars$debs, lower=0, upper=0.5,control=list(ftol=1e-15), quiet=TRUE))
		c(nls_result$value, nls_result$par[1], nls_result$par[2], nls_result$convergence)
	}
	stopImplicitCluster()
	colnames(nls_values) <- c("obj_fn","SS_parm","SD_parm","conv")
	nls_values <- data.frame(nls_values)

	best_init <- which.min(nls_values$obj_fn)
	nls_result <- nls_values[best_init,2:3]
	nls_coefs <- c(nls_result[1,1],nls_result[1,2])
	nls_best_init <- c(init_cond_grid[best_init,][1,1],init_cond_grid[best_init,][1,2])
	write.csv(nls_best_init,file="best_risk_negexp_nls_inits.csv")
}
if(find_best_nls_parms == 0) {
	nls_init <- read.csv("best_risk_negexp_nls_inits.csv")
	nls_result <- invisible(spg(par=c(nls_init[1,2],nls_init[2,2]), fn=objective, rate=risk, S=nls_risk_xvars$sats, D=nls_risk_xvars$debs, lower=0, upper=0.5,control=list(ftol=1e-15), quiet=TRUE))
	 nls_coefs <- c(nls_result$par[1],nls_result$par[2])
}

png(width=400,height=400,filename="../images/risk_fit_plot.png")
nls_fitplot(nls_risk_xvars,c(nls_coefs[1],nls_coefs[2]),series$year,risk,"Collision rate equation calibration","collision rate")
dev.off()

nls_coefs <- data.frame(S2=nls_coefs[1],SD=nls_coefs[2])

write.csv(nls_coefs,file="calibrated_risk_eqn_coefs.csv")

### linearized risk equation calibration
# lin_risk_xvars <- subset(series,select=c(payloads_in_orbit,debris))
# #lin_risk_xvars$log_sats <- log(lin_risk_xvars$payloads_in_orbit)
# #colnames(lin_risk_xvars) <- c("sats","debs","log_sats")
# colnames(lin_risk_xvars) <- c("sats","debs")
# linearized_risk <- log(lin_risk_xvars$sats - risk)
# linearized_risk <- log(1 - (risk/lin_risk_xvars$sats))

# test <- lm(linearized_risk ~ -1 + ., data=lin_risk_xvars)
# lin_risk_coefs <- coef(test)

# # plot fit
# lin_risk_xvars_mat <- as.matrix(lin_risk_xvars)
# fitplot(lin_risk_xvars_mat,lin_risk_coefs,series$year,linearized_risk,"Linearized collision rate")

# #m1 <- glmnet(x=lin_risk_xvars_mat,y=linearized_risk, alpha=0, lower.limits=c(-Inf, -Inf, 0), upper.limits=c(0, 0, Inf), lambda=cv.glmnet(x=lin_risk_xvars_mat,y=linearized_risk,alpha=0)$lambda.min,intercept=FALSE)
# m1 <- glmnet(x=lin_risk_xvars_mat,y=linearized_risk, alpha=0, lambda=cv.glmnet(x=lin_risk_xvars_mat,y=linearized_risk,alpha=0)$lambda.min,intercept=FALSE)
# coef(m1)

#####
# residual (semi-parametric) bootstrap for estimated nls parameters. grid searches same interval as original model estimate. 
#####
if(physics_bootstrap==1){
	B <- 1000 # number of bootstrap resamples

	nls_risk_resids <- risk - nls_risk_xvars[,1]*(1 - exp(-nls_coefs[1,1]*nls_risk_xvars[,1] -nls_coefs[1,2]*nls_risk_xvars[,2]))

	# generate bootstrap resamples of residuals 
	resample_order <- replicate(B,sample(c(1:length(nls_risk_resids)),length(nls_risk_resids),replace = TRUE))	
	nls_risk_resids_resamples <- matrix(nls_risk_resids[resample_order],ncol=B,byrow=FALSE)
	#replicate(B,sample(nls_risk_resids,length(nls_risk_resids),replace = TRUE))
	nls_risk_model_output <- nls_risk_xvars[,1]*(1 - exp(-nls_coefs[1,1]*nls_risk_xvars[,1] -nls_coefs[1,2]*nls_risk_xvars[,2]))
	nls_risk_bootstrap_dv <- nls_risk_model_output + nls_risk_resids_resamples

	nls_risk_bootstrap_coefs <- matrix(-1,nrow=B,ncol=2)

	registerDoParallel(cores=3)
	for(b in 1:B) {
		# search over a grid of initial points and select the best one
		init_risk_bootstrap_cond_element <- seq(from=1e-9,to=1e-6,length.out=10)
		init_risk_bootstrap_cond_grid <- expand.grid(init_risk_bootstrap_cond_element,init_risk_bootstrap_cond_element)
		nls_risk_bootstrap_values <- as.data.frame(matrix(1e+25,nrow=nrow(init_risk_bootstrap_cond_grid),ncol=4))
		colnames(nls_risk_bootstrap_values) <- c("obj_fn","SS_parm","SD_parm","conv")

		nls_risk_bootstrap_values <- foreach(i=1:nrow(init_risk_bootstrap_cond_grid), .export=ls(), .inorder=TRUE, .combine=rbind) %dopar% {
			nls_risk_bootstrap_result <- invisible(spg(par=c(init_risk_bootstrap_cond_grid[i,1],init_risk_bootstrap_cond_grid[i,2]), fn=objective, rate=nls_risk_bootstrap_dv[,b], S=nls_risk_xvars$sats, D=nls_risk_xvars$debs, lower=0, upper=0.5,control=list(ftol=1e-15), quiet=TRUE))
			c(nls_risk_bootstrap_result$value, nls_risk_bootstrap_result$par[1], nls_risk_bootstrap_result$par[2], nls_risk_bootstrap_result$convergence)
		}
		colnames(nls_risk_bootstrap_values) <- c("obj_fn","SS_parm","SD_parm","conv")
		nls_risk_bootstrap_values <- data.frame(nls_risk_bootstrap_values)

		best_init <- which.min(nls_risk_bootstrap_values$obj_fn)
		nls_risk_bootstrap_result <- nls_risk_bootstrap_values[best_init,2:3]
		nls_risk_bootstrap_coefs[b,] <- c(nls_risk_bootstrap_result[1,1],nls_risk_bootstrap_result[1,2])
	}
	stopImplicitCluster()

	nls_risk_bootstrap_coefs <- data.frame(nls_risk_bootstrap_coefs)
	colnames(nls_risk_bootstrap_coefs) <- c("S2","SD")

	# write out bootstrapped coefficients
	write.csv(nls_risk_bootstrap_coefs,file="bootstrapped_risk_eqn_coefs.csv")

	# plot bootstrapped coefficient distribution

	S2_hist <- ggplot(data=nls_risk_bootstrap_coefs, aes(x=S2)) + 
				geom_histogram(aes(y=..density..), binwidth=2e-8, colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
				xlab("Satellite-satellite collision parameter") +
				geom_vline(xintercept=nls_coefs[1,1], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(nls_risk_bootstrap_coefs$S2), linetype="dashed") +
				ggtitle(paste0("Distribution of bootstrapped satellite collision parameters: ",B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
				theme_bw()
	SD_hist <- ggplot(data=nls_risk_bootstrap_coefs, aes(x=SD)) + 
				geom_histogram(aes(y=..density..), binwidth=5e-9, colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
				xlab("Satellite-debris collision parameter") +
				geom_vline(xintercept=nls_coefs[1,2], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(nls_risk_bootstrap_coefs$SD), linetype="dashed") +
				ggtitle("\n") +
				theme_bw()
	grid.arrange(S2_hist,SD_hist,ncol=2)

	png(width=600,height=600,filename="../images/collision_rate_parameter_bootstrap_plot.png")
	grid.arrange(S2_hist,SD_hist,ncol=2)
	dev.off()

	nls_risk_bootstrap_dv_dfrm <- data.frame(nls_risk_bootstrap_dv)
	colnames(nls_risk_bootstrap_dv_dfrm) <- paste0("bootstrapped_dv_",seq(from=1,to=B,by=1))

	nls_risk_fitline <- nls_risk_model_output

	nls_risk_fit <- data.frame(year=series$year, model_fit=nls_risk_fitline, nls_risk_bootstrap_dv)
	nls_risk_fit_long <- reshape(nls_risk_fit, idvar="year", times=(colnames(nls_risk_fit)[-1]), direction="long", v.names="bootstrap_output", varying=(colnames(nls_risk_fit)[-1]))
	nls_risk_fit_long$model_fit <- nls_risk_fit$model_fit
	nls_risk_fit_long$truth <- risk
	rownames(nls_risk_fit_long) <- NULL

	nls_risk_bootstrap_plot_base <- ggplot(data=nls_risk_fit_long, aes(x=year))
	nls_risk_bootstrap_plot <- nls_risk_bootstrap_plot_base + geom_line(aes(y=bootstrap_output,group=as.factor(time)), size=0.5, alpha=0.2, color="light blue") + 
							geom_line(aes(y=model_fit), color="blue") +
							geom_line(aes(y=truth), size=1) +
							ylab("Expected number of collisions") +
							theme_minimal() + ggtitle(paste0("Residual bootstrap simulations of collision rate model: ",B," draws\n(blue: model fit, light blue: bootstrap fits, black: truth)"))
	nls_risk_bootstrap_plot


	png(width=450,height=400,filename="../images/collision_rate_bootstrap_prediction_plot.png")
	nls_risk_bootstrap_plot
	dev.off()
}

### create more variables for debris modeling
Sfrac <- 1 - exp(-nls_coefs[1,1]*dfrm$payloads_in_orbit)
Dfrac <- 1 - exp(-nls_coefs[1,2]*dfrm$debris)
dfrm$SSfrags <- Sfrac*dfrm$payloads_in_orbit
dfrm$SDfrags <- Dfrac*dfrm$payloads_in_orbit

series <- cbind(dfrm[-(dim(dfrm)[1]),],D_next,S_next,risk)
dfit_mat <- as.matrix(subset(series,select=c(payloads_in_orbit,payloads_decayed,debris,launch_successes,launch_failures,num_destr_asat,SSfrags,SDfrags) ) )
dfit_mat_elnet <- as.matrix(subset(dfit_mat,select=c(debris,launch_successes,num_destr_asat,SSfrags,SDfrags) ) )

### debris equation calibration
# Constrained ridge fit
m2 <- glmnet(x=dfit_mat_elnet,y=D_next, alpha=0, lower.limits=c(0, 0, 0, 0, 0), upper.limits=c(1, Inf, Inf, Inf, Inf), lambda=cv.glmnet(x=dfit_mat_elnet,y=D_next,alpha=0)$lambda.min,intercept=FALSE)
coef(m2)

# plot fit
m2xvars <- as.matrix(cbind(subset(series,select=colnames(dfit_mat_elnet))))
m2coefs <- coef(m2)[-1]
names(m2coefs) <- colnames(dfit_mat_elnet)
fitplot(m2xvars,m2coefs,series$year,D_next,"Debris law of motion calibration","debris")
dev.off()

# plot fit
png(width=400,height=400,filename="../images/debris_lom_fit_plot.png")
fitplot(m2xvars,m2coefs,series$year,D_next,"Debris law of motion calibration","debris")
dev.off()

write.csv(m2coefs,file="calibrated_debris_lom_coefs.csv")

#####
# residual (semi-parametric) bootstrap for estimated debris model parameters
#####

if(physics_bootstrap==1) {
	ridge_debris_resids <- D_next - m2xvars%*%m2coefs

	# generate bootstrap resamples of residuals 
	ridge_debris_resids_resamples <-  matrix(ridge_debris_resids[resample_order],ncol=B,byrow=FALSE)
	#replicate(B,sample(ridge_debris_resids,length(ridge_debris_resids),replace = TRUE))
	ridge_debris_model_output <- m2xvars%*%m2coefs
	ridge_debris_bootstrap_dv <- ridge_debris_model_output%*%matrix(1,ncol=B,nrow=1) + ridge_debris_resids_resamples

	ridge_debris_bootstrap_coefs <- matrix(-1,nrow=B,ncol=ncol(m2xvars))

	registerDoParallel(cores=3)
	ridge_debris_bootstrap_coefs <- foreach(b=1:B, .export=ls(), .inorder=TRUE, .combine=rbind) %dopar% {
		bootstrap_dfrm <- dfrm
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
	write.csv(ridge_debris_bootstrap_coefs,file="bootstrapped_debris_lom_coefs.csv")

	# draw histogram of coefficient distributions.
	title_hist <- ggtitle("\n")
	d_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=(1-debris))) + 
				geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
				xlab("Debris decay") +
				geom_vline(xintercept=(1-m2coefs[1]), linetype="dashed", color="blue") +
				geom_vline(xintercept=(1-mean(ridge_debris_bootstrap_coefs$debris)), linetype="dashed") +
				ggtitle(paste0("Distribution of bootstrapped debris law of motion parameters: ", B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
				theme_bw()
	m_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=launch_successes)) + 
				geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
				xlab("Launch debris") +
				geom_vline(xintercept=m2coefs[2], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$launch_successes), linetype="dashed") +
				ggtitle("\n") +
				theme_bw()
	gamma_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=num_destr_asat)) + 
				geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
				xlab(paste0("Fragments from \nanti-satellite missile tests")) +
				geom_vline(xintercept=m2coefs[3], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$num_destr_asat), linetype="dashed") +
				ggtitle("\n") +
				theme_bw()
	beta_SS_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SSfrags)) + 
				geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
				xlab(paste0("Fragments from \nsatellite-satellite collisions")) +
				geom_vline(xintercept=m2coefs[4], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SSfrags), linetype="dashed") +
				ggtitle("\n") +
				theme_bw()
	beta_SD_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SDfrags)) + 
				geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
				geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
				xlab(paste0("Fragments from \nsatellite-debris collisions")) +
				geom_vline(xintercept=m2coefs[5], linetype="dashed", color="blue") +
				geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SDfrags), linetype="dashed") +
				ggtitle("\n") +
				theme_bw()

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

	ridge_debris_bootstrap_plot_base <- ggplot(data=ridge_debris_fit_long, aes(x=year))
	ridge_debris_bootstrap_plot <- ridge_debris_bootstrap_plot_base + geom_line(aes(y=bootstrap_output,group=as.factor(time)), size=0.5, alpha=0.2, color="light blue") + 
							ylab("Debris in LEO") +
							geom_line(aes(y=model_fit), color="blue") +
							geom_line(aes(y=truth), size=1) +
							theme_minimal() + ggtitle(paste0("Residual bootstrap simulations of debris model: ",B," draws\n(blue: model fit, light blue: bootstrap fits, black: truth)"))
	ridge_debris_bootstrap_plot

	png(width=400,height=400,filename="../images/debris_bootstrap_prediction_plot.png")
	ridge_debris_bootstrap_plot
	dev.off()
}

