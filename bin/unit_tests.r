##### Unit tests for tihs project

#####
# Testing the optimal policy solver
#####
rm(list=ls())
library(data.table)
library(rootSolve)
library(gridExtra)
library(ggplot2)
library(viridis)
library(doParallel)
library(progress)
library(plot3D)
library(reshape2)
library(fields)
source("simulation_functions.r")
source("equations.r")
source("simulation_algorithms.r")

# helper function to calculate residuals and plot errors
fitplot <- function(betas,xvars,yvar,title) {
	fitline <- as.vector(xvars%*%betas)
	colnames(fitline)=NULL
	index <- c(1:length(fitline))
	fit <- data.frame(index=index,fit=fitline,truth=yvar,error=(fitline-yvar))

	plot_base <- ggplot(data=fit, aes(x=index))
	plot_fitplot <- plot_base + geom_line(aes(y=truth),size=1.1) +
							geom_line(aes(y=fit),size=0.9,linetype="dashed", color="blue") +
							theme_minimal() + ggtitle(paste(title))
	plot_error <- plot_base + geom_line(aes(y=error),size=0.9) +
						geom_hline(yintercept=0,linetype="dashed") +
						theme_minimal()

	grid.arrange(plot_fitplot,plot_error,nrow=2)
}

# build grid, generate guesses, initialize dynamic_vfi_solver output list
gridsize <- 64
gridlist <- build_grid(gridmin=0, gridmax=17500, gridsize,1)
vguess <- matrix(gridlist$igrid$sats,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
dvs_output <- list()

# define T, sequence of p, F, and asats
T <- 10
p <- rep(1,length=T)
F <- seq(from=15,to=13,length=T)
asats <- rep(0,length=T)

# define physical parameters and discount rate+factor
risk_cal <- read.csv("../data/calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("../data/calibrated_debris_lom_coefs.csv")
observed_time_series <- read.csv("../data/ST_ESA_series.csv")

start_year <- 2005
S0 <- observed_time_series$payloads_in_orbit[which(observed_time_series$year==start_year)]
D0 <- observed_time_series$debris[which(observed_time_series$year==start_year)]
avg_sat_decay <- mean(observed_time_series$payloads_decayed[which(observed_time_series$year>=start_year)]/observed_time_series$payloads_in_orbit[which(observed_time_series$year>=start_year)])

aS <- risk_cal[1,2]
aD <- risk_cal[2,2]
aSS <- risk_cal[3,2]
aSD <- risk_cal[4,2]
aDD <- risk_cal[5,2]

aDDbDD <- deblom_cal[7,2]
bSS <- deblom_cal[5,2]
bSD <- deblom_cal[6,2]
d <- deblom_cal[2,2]
Z_coef <- deblom_cal[1,2]
m <- deblom_cal[3,2]
asat_coef <- deblom_cal[4,2]
#asats <- observed_time_series$num_destr_asat[which(observed_time_series$year>=start_year)]

discount_rate <- 0.05#0.6660755
discount_fac <- 1/(1+discount_rate)

# Check that objective function makes sense
k=50
fleet_preval(X=0,S=gridpanel$S[k],D=gridpanel$D[k],value_fn=gridpanel$V,asats=asats,t=T,p=p,F=F,igrid=gridlist$igrid)

# Check that path solver works in all periods
registerDoParallel(cores=4)
for(i in T:1){
	dvs_output[[i]] <- dynamic_vfi_solver(gridpanel,igrid=gridlist$igrid,asats,i,T,p,F)
	vguess <- matrix(dvs_output[[i]]$optimal_fleet_vfn,nrow=gridsize,ncol=gridsize)
	lpguess <- matrix(dvs_output[[i]]$optimal_launch_pfn,nrow=gridsize,ncol=gridsize)
	gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
}
stopImplicitCluster()

# bind the list of solved policies into a long dataframe
policy_path <- rbindlist(dvs_output)
optimal_launch_pfn <- as.vector(policy_path$optimal_launch_pfn)

# Thin plate spline to fit the policy function
tps_x <- as.matrix(cbind(policy_path$satellites,policy_path$debris,policy_path$F))
tps_y <- as.matrix(optimal_launch_pfn)
tps_model <- Tps(x=tps_x,Y=tps_y)
surface(tps_model)

predict(tps_model,x=cbind(5.1e+3,3.9e+7,F[5]))

# compare tps time series to solved time series

### function to begin an optimal launch sequence at a given time
simulate_optimal_path <- function(p,F,discount_rate,T,...) {
	fe_eqm <- p/F - discount_rate
	asats_inf <- rep(0,length=T)
	launch_constraint_inf <- rep(1e+10,length=T)
	opt_path <- fp_tsgen(0,0,T,fe_eqm,launch_constraint_inf,asats_inf,p,F)

	return(opt_path)
}

tps_opt_path <- function(S0,D0,p,F,tps_model,asats_seq) {
	times <- seq(from=1,to=T,by=1)	
	sat_seq <- rep(0,length=T)
	deb_seq <- rep(0,length=T)
	profit_seq <- rep(0,length=T)
	X <- rep(-1,length=T)

	sat_seq[1] <- S0
	deb_seq[1] <- D0
	X[1] <- predict(tps_model,x=cbind(sat_seq[1],deb_seq[1],F[1]))
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1,p,F)

	for(k in 2:T) {
		sat_seq[k] <- S_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)])
		deb_seq[k] <- D_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)],asats_seq[(k-1)])
		X[k] <- predict(tps_model,x=cbind(sat_seq[k],deb_seq[k],F[k]))
		X[k] <- ifelse(X[k]<0,0,X[k])
		profit_seq[k] <- one_p_return(X[k],sat_seq[k],k,p,F)*(discount_fac^(times[(k-1)]))
	}
	deb_seq[is.na(deb_seq)] <- max(!is.na(deb_seq))
	profit_seq[T] <- fleet_ssval_T(X[T],sat_seq[T],T,p,F)
	losses <- L(sat_seq,deb_seq)
	values <- as.data.frame(cbind(times,X,sat_seq,deb_seq,profit_seq,losses))
	colnames(values) <- c("time","launches","satellites","debris","fleet_pv","collision_rate")
	return(values)
}

opt_path <- simulate_optimal_path(p,F,discount_rate,T)
tps_path <- tps_opt_path(0,0,p,F,tps_model,asats)
View(opt_path)
View(tps_path)
