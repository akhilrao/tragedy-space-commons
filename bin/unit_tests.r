##### Unit tests for this project

#####
# Testing the optimal policy solver
#####
rm(list=ls())
library(pracma)
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
system(sprintf("taskset -p 0xffffffff %d", Sys.getpid())) # Adjusts the R session's affinity mask from 1 to f, allowing the R process to use all cores.

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
gridsize <- 32
gridlist <- build_grid(gridmin=0, gridmax=35000, gridsize, cheby=1)
### shift grid down to zero if chebysheving it moved it up - should be unnecessary with expanded chebysehv array (secant factor)
if(min(gridlist$base_piece)>0) {
	gridlist$base_piece <- gridlist$base_piece - min(gridlist$base_piece)
	gridlist$igrid <- gridlist$igrid - min(gridlist$igrid)
}
vguess <- matrix(gridlist$igrid$sats,nrow=gridsize,ncol=gridsize)
lpguess <- matrix(0,nrow=gridsize,ncol=gridsize)
gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
dvs_output <- list()

# define T, sequence of p, F, and asats
T <- 100
p <- rep(1,length=T)
F <- c(rep(13,length=T/2),rep(10,length=T/2))
#F <- seq(from=15,to=15,length.out=T)
#F <- rep(13,length=T)
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

# Check that path solver works in all periods
total.grid.time <- proc.time()[3]
registerDoParallel(cores=32)
for(i in T:1){
	dvs_output[[i]] <- dynamic_vfi_solver(gridpanel,igrid=gridlist$igrid,asats,i,T,p,F)
	vguess <- matrix(dvs_output[[i]]$optimal_fleet_vfn,nrow=gridsize,ncol=gridsize)
	lpguess <- matrix(dvs_output[[i]]$optimal_launch_pfn,nrow=gridsize,ncol=gridsize)
	gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
	dev.off()
}
stopImplicitCluster()
cat(paste0("\n Done. Total grid compute time taken: ",round(proc.time()[3] - total.grid.time,3)," seconds"))

# bind the list of solved policies into a long dataframe
policy_path <- rbindlist(dvs_output)

# compare tps time series to solved time series
### function to begin an optimal launch sequence at a given time
simulate_optimal_path <- function(p,F,discount_rate,T,...) {
	fe_eqm <- p/F - discount_rate
	asats_inf <- rep(0,length=T)
	launch_constraint_inf <- rep(1e+10,length=T)
	opt_path <- fp_tsgen(0,0,T,fe_eqm,launch_constraint_inf,asats_inf,p,F)

	return(opt_path)
}

tps_opt_path <- function(S0,D0,p,F,policy_path,asats_seq,igrid,ncores) {
	times <- seq(from=1,to=T,by=1)	
	sat_seq <- rep(0,length=T)
	deb_seq <- rep(0,length=T)
	profit_seq <- rep(0,length=T)
	X <- rep(-1,length=T)

	sat_seq[1] <- S0
	deb_seq[1] <- D0
	optimal_launch_pfn <- as.vector(policy_path$optimal_launch_pfn)

	## Thin plate splines to fit the policy functions
	# current_cost <- which(igrid$F==F[1])
	# tps_x <- as.matrix(cbind(policy_path$satellites[current_cost],policy_path$debris[current_cost]))
	# tps_y <- as.matrix(optimal_launch_pfn[current_cost])
	# cat(paste0("\nEstimating spline interpolant of period 1 policy function..."))
	# tps_model <- Tps(x=tps_x,Y=tps_y)
	# cat(paste0("\n Done."))

	spline_list <- list()

	s.tm <- proc.time()[3]
	cat(paste0("\nEstimating spline interpolants of policy functions..."))
	spline_list <- foreach(k=1:T, .export=ls(), .inorder=TRUE) %dopar% {
			# cat(paste0("\nEstimating spline interpolant of period ", k, " policy function..."))
			current_cost <- which(igrid$F==F[k])
			tps_x <- as.matrix(cbind(policy_path$satellites[current_cost],policy_path$debris[current_cost]))
			tps_y <- as.matrix(optimal_launch_pfn[current_cost])
			tps_model <- Tps(x=tps_x,Y=tps_y)
			return(tps_model)
		}
	cat(paste0("\n Done. Total time taken: ",round(proc.time()[3] - s.tm,3)," seconds"))

	s.tm <- proc.time()[3]
	cat(paste0("\nGenerating policy time path..."))
	# s.tm <- proc.time()[3]
	# cat(paste0("\nGenerating period 1 policy..."))
	X[1] <- predict(spline_list[[1]],x=cbind(sat_seq[1],deb_seq[1]))
	X[1] <- ifelse(X[1]<0,0,X[1])
	# cat(paste0("\n Done. Time taken: ",round(proc.time()[3] - s.tm,3)))
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1,p,F)

	for(k in 2:T) {
		# # re-estimate spline for current period
		# current_cost <- which(igrid$F==F[k])
		# tps_x <- as.matrix(cbind(policy_path$satellites[current_cost],policy_path$debris[current_cost]))
		# tps_y <- as.matrix(optimal_launch_pfn[current_cost])
		# s.tm <- proc.time()[3]
		# cat(paste0("\nEstimating spline interpolant of period ", k, " policy function..."))
		# tps_model <- Tps(x=tps_x,Y=tps_y)
		# cat(paste0("\n Done. Time taken: ",round(proc.time()[3] - s.tm,3)))
		# compute next state and interpolated policy
		sat_seq[k] <- S_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)])
		deb_seq[k] <- D_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)],asats_seq[(k-1)])
		# s.tm <- proc.time()[3]
		# cat(paste0("\nGenerating period ", k," policy..."))
		X[k] <- predict(spline_list[[k]],x=cbind(sat_seq[k],deb_seq[k]))
		X[k] <- ifelse(X[k]<0,0,X[k])
		# cat(paste0("\n Done. Time taken: ",round(proc.time()[3] - s.tm,3)))
		profit_seq[k] <- one_p_return(X[k],sat_seq[k],k,p,F)*(discount_fac^(times[(k-1)]))
	}
	cat(paste0("\n Done. Total time taken: ",round(proc.time()[3] - s.tm,3)," seconds"))
	deb_seq[is.na(deb_seq)] <- max(!is.na(deb_seq))
	profit_seq[T] <- fleet_ssval_T(X[T],sat_seq[T],T,p,F)
	losses <- L(sat_seq,deb_seq)
	values <- as.data.frame(cbind(times,X,sat_seq,deb_seq,profit_seq,losses,p,F))
	colnames(values) <- c("time","launches","satellites","debris","fleet_pv","collision_rate","returns","costs")
	return(values)
}

linint_opt_path <- function(S0,D0,p,F,policy_lookup,asats_seq,igrid) {
	times <- seq(from=1,to=T,by=1)	
	sat_seq <- rep(0,length=T)
	deb_seq <- rep(0,length=T)
	profit_seq <- rep(0,length=T)
	X <- rep(-1,length=T)

	sat_seq[1] <- S0
	deb_seq[1] <- D0
	next_state <- c(sat_seq[1],deb_seq[1])
	current_cost <- which(igrid$F==F[1])
	X[1] <- interpolate(next_state,igrid[current_cost,],policy_lookup[current_cost])
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1,p,F)

	for(k in 2:T) {
		sat_seq[k] <- S_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)])
		deb_seq[k] <- D_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)],asats_seq[(k-1)])	
		next_state <- c(sat_seq[k],deb_seq[k])
		current_cost <- which(igrid$F==F[k])
		X[k] <- interpolate(next_state,igrid[current_cost,],policy_lookup[current_cost])
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

grid_lookup <- data.frame(sats=policy_path$satellites,debs=policy_path$debris,F=policy_path$F)

tps_path <- tps_opt_path(0,0,p,F,policy_path,asats,grid_lookup,4)
View(tps_path)
li_path <- linint_opt_path(0,0,p,F,policy_path$optimal_launch_pfn,asats,grid_lookup)
View(li_path)

tps_path_base <- ggplot(tps_path, aes(x=time)) + xlab("Time") + theme_minimal()
tps_launches <- tps_path_base + geom_line(aes(y=launches), size=1) + ylab("Launches")
tps_sats <- tps_path_base + geom_line(aes(y=satellites), size=1) + ylab("Satellites")
tps_debs <- tps_path_base + geom_line(aes(y=debris), size=1) + ylab("Debris")
tps_risk <- tps_path_base + geom_line(aes(y=collision_rate), size=1) + ylab("Collision risk")
tps_returns <- tps_path_base + geom_line(aes(y=returns), size=1) + ylab("Returns")
tps_costs <- tps_path_base + geom_line(aes(y=costs), size=1) + ylab("Costs")

grid.arrange(tps_launches,tps_sats,tps_debs,tps_risk,tps_returns,tps_costs,nrow=3,ncol=2)

png(file=paste0("timepath_",gridsize,"_pt_basegrid.png"))
grid.arrange(tps_launches,tps_sats,tps_debs,tps_risk,tps_returns,tps_costs,nrow=3,ncol=2)
dev.off()
