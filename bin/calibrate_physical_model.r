#####################################################################
# Script to take a series of orbital stock levels and calibrate a physical model of debris evolution and collision risk
#####################################################################

setwd("../data/")
set.seed(501) # for bootstrap samples

dfrm <- read.csv("ST_stock_series.csv")
# limit to years after 1990
#dfrm <- dfrm[which(dfrm$year>=1990),]
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]

##### Create additional RHS variables for calibration regressions
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
# fitplot(sat_xvars,sat_coefs,series$year,S_next,"Satellite accounting model","satellites")
# dev.off()

### Satellite equation calibration. this is the statistical equation: calculating decay rate implied by including collisions. Produces a positive, large, and insignificant coefficient on risk, in addition to about 0.5 and significant coefficient on launch_successes. Implied satellite time on-orbit is close to 30 years, though, which seems reasonable.
# Approach: Fit model to estimate decay coefficient, then set other coefficients to (-decay_coefficient, 1). The fit from this is not much worse than the fit with the estimated coefficients.
sfit_dfrm <- subset(series,select=c(payloads_in_orbit,risk,launch_successes) )

sat_ols <- lm(S_next ~ -1 + payloads_in_orbit + risk + launch_successes, data=sfit_dfrm)
summary(sat_ols)
sat_parms <- coef(sat_ols)

# plot fit
sat_xvars <- as.matrix(cbind(subset(sfit_dfrm,select=colnames(sfit_dfrm))))
sat_coefs <- c(sat_parms[1],-1*sat_parms[1],1)

write.csv(sat_parms[1],file="calibrated_satellite_lom_coefs.csv")

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

nls_coefs <- data.frame(S2=nls_coefs[1],SD=nls_coefs[2])

write.csv(nls_coefs,file="calibrated_risk_eqn_coefs.csv")


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

write.csv(m2coefs,file="calibrated_debris_lom_coefs.csv")

setwd("../bin/") # set working directory back to bin for the next stage of the main analysis

if(physics_bootstrap==1){
	source("physics_bootstrap.r") # this regenerates Extended Data figure 7
}


##### Generate figures for Extended Data section

# Extended Data figure 2

ed_fig2_a <- fitplot_noerror(sat_xvars,sat_coefs,series$year,S_next,"Satellite law of motion calibration","satellites")
ed_fig2_b <- fitplot_noerror(m2xvars,m2coefs,series$year,D_next,"Debris law of motion calibration","debris")
ed_fig2_c <- nls_fitplot_noerror(nls_risk_xvars,c(nls_coefs[1,1],nls_coefs[1,2]),series$year,risk,"Collision rate physical equation calibration","collision rate")

png(width=950,height=400,filename="../images/extended_data_figure_2.png")
plot_grid(ed_fig2_a, ed_fig2_b, ed_fig2_c, labels=c("a", "b", "c"), align="h", ncol=3, label_size=20)
dev.off()
