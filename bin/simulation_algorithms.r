##### algorithms to be used in deterministic satellite-debris model dynamic programming simulation script

# fleet planner time series generator
fp_tsgen <- function(S,D,T,fe_eqm,launch_con,asats_seq,p,F,...) {
	ctm <- proc.time()[3]
	parms <- c(rep(15,length=T),S,D)
	launch_path <- optim(par=parms[1:T],fn=fhvf, S=parms[(T+1)], D=parms[(T+2)], asats_seq=asats_seq, launch_con=launch_con, T=T,p=p,F=F, control=list(fnscale=-1, pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1e+14)$par
	launch_path_clock <- proc.time()[3] - ctm

	ctm <- proc.time()[3]
	time_series <- seriesgen_ts(launch_path,S,D,T,asats_seq,p,F)
	propagation_clock <- proc.time()[3] - ctm

	print("Fleet planner time series generated.") 
	print(paste0("Time to compute launch path: ", launch_path_clock, " seconds."))
	print(paste0("Time to propagate stocks: ", propagation_clock, " seconds."))

	return(time_series)
}

# open access time series generator
oa_tsgen <- function(S,D,T,fe_eqm,launch_con,asats,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	profit_seq <- rep(0,length=(T+1))
	discounted_profit_seq <- rep(0,length=(T+1))
	for(t in 1:T ) {
		X <- uniroot.all(eqmcond,c(0,1e+6),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm[t+1],asats=asats[t])
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		if(launch_path[t]>launch_con[t]) {launch_path[t] <- launch_con[t]}
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t],asats[t])
		profit_seq[t] <- one_p_return(launch_path[t],sat_path[t],t,p,F)
		discounted_profit_seq[t] <- profit_seq[t]*discount_fac^(t-1)
	}
	clock <- proc.time()[3] - ctm
	#launch_path[T+1] <- uniroot.all(eqmcond,c(0,1e+6),S=sat_path[T+1],D=deb_path[T+1],fe_eqm=fe_eqm[T+1],asats=asats[T+1])
	print(uniroot.all(eqmcond,c(0,1e+6),S=sat_path[T+1],D=deb_path[T+1],fe_eqm=fe_eqm[T+1],asats=asats[T+1]))

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_flowv=profit_seq,fleet_pv=discounted_profit_seq,collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

### open access policy algorithm - full serial (still fast - just rootfinding with no integration)
oapolicy <- function(igrid,fe_eqm,t,asats,p,F,...) {
ptm <- proc.time()
	launches <- rep(0,length=dim(igrid)[1])
	fe_eqm_vec <- rep(fe_eqm,length=dim(igrid)[1])
	p_vec <- rep(p,length=dim(igrid)[1])
	F_vec <- rep(F,length=dim(igrid)[1])
	for(j in 1:dim(igrid)[1]) {
		X <- uniroot.all(eqmcond,c(0,1e+9),S=igrid$sats[j],D=(igrid$debs[j]),fe_eqm=fe_eqm[t+1],asats=asats[t])
		launches[j] <- ifelse(length(X)==0,0,X)
	}

	launches[which(launches=="Inf")] <- -1
	fleet_size <- S_(launches,igrid$sats,igrid$debs)

	loss <- L(S_(launches,igrid$sats,igrid$debs),D_(launches,igrid$sats,igrid$debs,asats))
	results <- as.data.frame(cbind(igrid,launches,fleet_size,loss,fe_eqm_vec,p_vec,F_vec))
	colnames(results) <- c("satellites","debris","oa_launch_pfn","oa_fleet_size","loss_next","fe_eqm","p","F")
	print(paste("Total time taken for open access policy: ",round(((proc.time() - ptm)[3])/60,4)," minutes",sep=""))
	return(results)
}

### algorithm to compute open access value functions along a given returns and cost path
oa_pvfn_path_solver <- function(dvs_output,gridpanel,gridlist,asats,T,p,F,fe_eqm,ncores,...) {
	total.grid.time <- proc.time()[3]
	registerDoParallel(cores=ncores)
	base_grid <- unique(gridlist$igrid[,1])
	fe_eqm_padded <- c(fe_eqm,fe_eqm[T])
	for(i in T:1){
		dvs_output[[i]] <- oapolicy(igrid=gridlist$igrid,fe_eqm_padded,t=i,asats=asats,p=p[i],F=F[i])
		if(i==T) {
			tps_x <- as.matrix(cbind(gridlist$igrid$sats,gridlist$igrid$debs))
			tps_y <- as.matrix(gridpanel$V)
			tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y, lambda=0))
			spline_vfn_int <- as.vector(predict(tps_model,x=tps_x))
			spline_vfn_int_mat <- matrix(spline_vfn_int,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
			pfn <- matrix(dvs_output[[i]]$oa_launch_pfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
			dev.new(width=8,height=7,unit="in")
			par(mfrow=c(2,2))
			image2D(z=spline_vfn_int_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("value function interpolation"))
			image2D(z=spline_vfn_int_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("value function interpolation"))
			image2D(z=pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("policy function"))
			image2D(z=pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("policy function"))
			value_fn <- policy_eval_BI(igrid=gridlist$igrid,launch_policy=dvs_output[[i]]$oa_launch_pfn,value_fn=V_T,T=250,tps_model=tps_model,p_t=p[T],F_t=F[T],asats_t=0)
			dvs_output[[i]]$oa_fleet_vfn <- value_fn
		}
		if(i!=T) {
			tps_x <- as.matrix(cbind(gridlist$igrid$sats,gridlist$igrid$debs))
			tps_y <- as.matrix(dvs_output[[i+1]]$oa_fleet_vfn)
			spline_vfn_int <- as.vector(predict(tps_model,x=tps_x))
			spline_vfn_int_mat <- matrix(spline_vfn_int,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)			
			pfn <- matrix(dvs_output[[i]]$oa_launch_pfn,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
			dev.new(width=8,height=7,unit="in")
			par(mfrow=c(2,2))
			image2D(z=spline_vfn_int_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("value function interpolation"))
			image2D(z=spline_vfn_int_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("value function interpolation"))
			image2D(z=pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("policy function"))
			image2D(z=pfn,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=FALSE,main=c("policy function"))
			value_fn <- foreach(j=1:length(gridlist$igrid$sats), .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
				result <- suppressWarnings(fleet_preval_spline(X=dvs_output[[i]]$oa_launch_pfn[j],S=gridlist$igrid$sats[j],D=gridlist$igrid$debs[j],value_fn=value_fn,asats=asats,t=i,p=p,F=F,igrid=gridlist$igrid,tps_model=tps_model))
				result
			}
			dvs_output[[i]]$oa_fleet_vfn <- value_fn
		}
		dev.off()
	}
	stopImplicitCluster()
	cat(paste0("\n Done. Total grid compute time taken: ",round(proc.time()[3] - total.grid.time,3)," seconds"))
	return(dvs_output)
}

### Computes backwards induction from given terminal value with given policy to get value function
policy_eval_BI <- function(igrid,launch_policy,value_fn,T,tps_model,p_t,F_t,asats_t) {
	count <- 0
	p_vec <- rep(p_t,length=T)
	F_vec <- rep(F_t,length=T)
	asats_vec <- rep(0,length=T) 
	tot.time <- proc.time()[3]
	## make progress bar
	BIpb <- progress_bar$new(format="Doing backwards induction [:bar] :percent",total=(T+1))
	BIpb$tick(0)
	while(count<=T) {
		value_fn <- foreach(i=1:length(launch_policy), .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			result <- fleet_preval_spline(X=launch_policy[i],S=igrid$sats[i],D=igrid$debs[i],value_fn=value_fn,asats=asats_vec,t=T,p=p_vec,F=F_vec,igrid=igrid,tps_model=tps_model)
			result
		}
		count <- count + 1
		BIpb$tick()
	}
	BIpb$terminate()
	tot.time <- round(round(proc.time()[3] - tot.time,3)/60,3)
	cat(paste0("\nTotal time taken for backwards induction: ", tot.time, " minutes.\n Average time per period: ", round(tot.time/length(launch_policy),3), " minutes.\n"))
	return(value_fn)
}

### fleet planner VFI algorithm: solve for policy assuming steady state reached or assuming a perfect-foresight path to steady state
dynamic_vfi_solver <- function(panel,igrid,asats,t,T,p,F,...) {
	# initialize hyperparameters
	panrows <- nrow(panel)
	gridmin <- min(igrid)
	gridmax <- max(igrid)
	n_grid_points <- length(unique(igrid[,1]))^2
	base_grid <- unique(igrid[,1])
	#dev.new(width=12,height=5,unit="in")
	dev.new(width=8,height=7,unit="in")
	par(mfrow=c(2,2))

	# initialize output objects
	newX <- rep(-1,length=panrows)
	result <- cbind(newX,newX,new)
	# initialize epsilon-delta and count
	ifelse(t==T, epsilon <- max(n_grid_points*1e-5,1e-3), epsilon <- max(n_grid_points*2e-6,1e-3)) # tighter epsilon for value function convergence in final period, looser epsilon for policy function convergence in prior periods.
	#ifelse(t==T, epsilon <- 10, epsilon <- 1) # for testing
	ifelse(t==T,panel$X <- panel$X, panel$X <- panel$X) #rnorm(length(panel$X),mean=10,sd=1))
	delta_old <- 0
	delta <- epsilon + 10
	delta2 <- 25
	count <- 0
	cat(paste("\n Stopping criteria: distance < ", epsilon, "\n", sep=""))

	# solver loop
	while(delta > epsilon) {

		## plot pfn and vfn
		plot_pfn_vfn(panel$V,panel$X,base_grid,c("Value function","Policy function"))

		## create spline interpolation model
		cat(paste0("\nEstimating spline interpolant of value function..."))
		tps_x <- as.matrix(cbind(panel$S,panel$D))
		tps_y <- as.matrix(panel$V)
		tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y, lambda=0))
		cat(paste0("\n Done."))

		spline_vfn_int <- as.vector(predict(tps_model,x=tps_x))
		spline_vfn_int_mat <- matrix(spline_vfn_int,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
		policy_mat <- matrix(panel$X,nrow=length(base_grid),ncol=length(base_grid),byrow=TRUE)
		# TODO: replace this with call to plot_pfn_vfn
		image2D(z=spline_vfn_int_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("value function interpolation"))
		image2D(z=policy_mat,x=base_grid,y=base_grid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=c("policy function"))

		t.tm <- proc.time()
		## maximization step
		cat(paste0("\n Starting maximization step, grid size = ", n_grid_points,": \n","."))
		m.tm <- proc.time()
		result <- foreach(k=1:panrows, .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			ctm <- proc.time()[3]
			ulim <- gridmax
			solution <- optim(par=panel$X[k], fn=fleet_preval_spline, S=panel$S[k], D=panel$D[k], value_fn=panel$V, asats=asats, t=t, p=p, F=F, igrid=igrid, tps_model=tps_model, control=list(fnscale=-1), method="L-BFGS-B", lower=0, upper=ulim)
			launch_rate <- solution$par[1]
			vfn <- solution$value
			clock <- proc.time()[3] - ctm
			result <- c(launch_rate,vfn,clock)
			return(result)
		}

		maxim.time <- (proc.time()-m.tm)[3]
		avg_cell_time <- mean(result[,3])
		cat(paste0("Finished. \n Average cell time: ",round(avg_cell_time,3)," seconds.\nTotal grid time: ",round((avg_cell_time*n_grid_points),3)," seconds \n"))
		rownames(newX) <- NULL

		## save new policy and new value
		newX <- result[,1]
		newV <- result[,2]

		## calculate distance (|V-newV| or |X-newX|) and update policy or value  
		if(t==T) {
			policy_delta <- max(abs((panel$X-newX)))
			ifelse(count==0, newV <- policy_eval_BI(igrid,newX,newV,T=10,tps_model,p[T],F[T],asats[T]), newV <- newV)
			ifelse(policy_delta<1, newV <- policy_eval_BI(igrid,newX,newV,T=min(count+1,75),tps_model,p[T],F[T],asats[T]), newV <- newV)
			cat(paste("\n Policy delta is ", policy_delta, sep=""))
			delta <- max(abs((panel$V-newV)))
			panel$V <- newV
			panel$X <- newX
		}
		if(t!=T) {
			delta <- max(abs((panel$X-newX)))
			panel$X <- newX			
		}

		## update old delta and calculate percentage change in delta by midpoint method		
		delta2 <- abs( (delta-delta_old)/(0.5*(delta+delta_old)) )
		delta_old <- delta

		## update count
		count <- count+1
	
		## calculate total time
		tot.time <- (proc.time() - t.tm)[3]
		## out
		cat(paste("\n Finished period ", t," iteration ", count,", distance is ", delta, "\n Stopping criteria: distance < ", epsilon, "\n", sep=""))
	}

	if(t!=T) {panel$V <- newV}
	# png(file=paste0("Value_policy_t_",t,".png"))
	# plot_pfn_vfn(panel$V,panel$X,base_grid,c("Value function","Policy function"))
	# dev.off()
	output <- as.data.frame(cbind(satellites=panel$S,debris=panel$D,optimal_launch_pfn=panel$X,optimal_fleet_vfn=panel$V,optimal_fleet_size=S_(panel$X,panel$S,panel$D),t=t,p=p[t],F=F[t]))
	colnames(output)[4] <- "optimal_fleet_vfn"
	return(output)
}

### algorithm to compute optimal policy and value functions along a given returns and cost path
opt_pvfn_path_solver <- function(dvs_output,gridpanel,gridsize,gridlist,asats,T,p,F,ncores,...) {
	total.grid.time <- proc.time()[3]
	registerDoParallel(cores=ncores)
	for(i in T:1){
		dvs_output[[i]] <- dynamic_vfi_solver(gridpanel,igrid=gridlist$igrid,asats,i,T,p,F)
		vguess <- matrix(dvs_output[[i]]$optimal_fleet_vfn,nrow=gridsize,ncol=gridsize)
		lpguess <- matrix(dvs_output[[i]]$optimal_launch_pfn,nrow=gridsize,ncol=gridsize)
		gridpanel <- grid_to_panel(gridlist,lpguess,vguess)
		dev.off()
	}
	stopImplicitCluster()
	cat(paste0("\n Done. Total grid compute time taken: ",round(proc.time()[3] - total.grid.time,3)," seconds"))
	return(dvs_output)
}

### function to begin an optimal finite-horizong launch sequence at a given time
simulate_optimal_path <- function(p,F,discount_rate,T,...) {
	fe_eqm <- p/F - discount_rate
	asats_inf <- rep(0,length=T)
	launch_constraint_inf <- rep(1e+10,length=T)
	opt_path <- fp_tsgen(0,0,T,fe_eqm,launch_constraint_inf,asats_inf,p,F)

	return(opt_path)
}

### algorithm to compute a time path using thin plate spline interpolation of solved policy functions
# S0, D0 are the initial conditions. t0 is the initial time, which is used to generate paths starting at different times.
tps_path_gen <- function(S0,D0,t0,p,F,policy_path,asats_seq,launchcon_seq,igrid,ncores,OPT,linear_policy_interp) {
	times <- seq(from=1,to=(T-t0),by=1)	
	sat_seq <- rep(0,length=(T-t0))
	deb_seq <- rep(0,length=(T-t0))
	profit_seq <- rep(0,length=(T-t0))
	discounted_profit_seq <- rep(0,length=(T-t0))
	fleet_npv_path <- rep(0,length=(T-t0))
	X <- rep(-1,length=(T-t0))

	if(length(launchcon_seq)==-1) {launchcon_seq <- rep(upper,length=(T-t0))} # -1 is a flag to set the constraint large enough that it never binds

	sat_seq[1] <- S0
	deb_seq[1] <- D0
	ifelse(OPT==1, launch_pfn <- as.vector(policy_path$optimal_launch_pfn), launch_pfn <- as.vector(policy_path$oa_launch_pfn))
	ifelse(OPT==1, fleet_vfn <- as.vector(policy_path$optimal_fleet_vfn), fleet_vfn <- as.vector(policy_path$oa_fleet_vfn))

	# Thin plate splines to fit the policy functions
	spline_list <- list()

	s.tm <- proc.time()[3]
	cat(paste0("\nEstimating spline interpolants of policy functions..."))
	spline_list <- foreach(k=1:(T-t0), .export=ls(), .inorder=TRUE) %dopar% {
			current_cost <- which(igrid$F==F[k])
			tps_x <- as.matrix(cbind(policy_path$satellites[current_cost],policy_path$debris[current_cost]))
			tps_y <- as.matrix(launch_pfn[current_cost])
			ifelse(linear_policy_interp==1,tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y,lambda=0)),tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y)))
			return(tps_model)
		}
	cat(paste0("\n Done. Total time taken: ",round(proc.time()[3] - s.tm,3)," seconds"))

	cat(paste0("\nEstimating spline interpolants of value functions..."))
	vfn_spline_list <- foreach(k=1:(T-t0), .export=ls(), .inorder=TRUE) %dopar% {
			current_cost <- which(igrid$F==F[k])
			tps_x <- as.matrix(cbind(policy_path$satellites[current_cost],policy_path$debris[current_cost]))
			tps_y <- as.matrix(fleet_vfn[current_cost])
			vfn_tps_model <- suppressWarnings(Tps(x=tps_x,Y=tps_y,lambda=0))
			return(vfn_tps_model)
		}
	cat(paste0("\n Done. Total time taken: ",round(proc.time()[3] - s.tm,3)," seconds"))

	s.tm <- proc.time()[3]
	cat(paste0("\nGenerating policy time path..."))
	X[1] <- predict(spline_list[[1]],x=cbind(sat_seq[1],deb_seq[1]))
	X[1] <- ifelse(X[1]<0,0,X[1])
	X[1] <- ifelse(X[1]>launchcon_seq[1],X[1]<-launchcon_seq[1],X[1]<-X[1])
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1,p,F)
	discounted_profit_seq[1] <- one_p_return(X[1],sat_seq[1],1,p,F)
	fleet_npv_path[1] <- predict(vfn_spline_list[[1]],x=cbind(sat_seq[1],deb_seq[1]))

	for(k in 2:(T-t0)) {
		current_clock_time <- t0 + k # need to calculate what the time is in the outside world for asats and launch constraint
		sat_seq[k] <- S_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)])
		deb_seq[k] <- D_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)],asats_seq[(current_clock_time-1)])
		X[k] <- predict(spline_list[[k]],x=cbind(sat_seq[k],deb_seq[k]))
		X[k] <- ifelse(X[k]<0,0,X[k])
		X[k] <- ifelse(X[k]>launchcon_seq[current_clock_time],X[k]<-launchcon_seq[current_clock_time],X[k]<-X[k])
		profit_seq[k] <- one_p_return(X[k],sat_seq[k],k,p,F)
		discounted_profit_seq[k] <- profit_seq[k]*(discount_fac^(times[(k-1)]))
		fleet_npv_path[k] <- predict(vfn_spline_list[[k]],x=cbind(sat_seq[k],deb_seq[k]))
	}
	cat(paste0("\n Done. Total time taken: ",round(proc.time()[3] - s.tm,3)," seconds"))
	deb_seq[is.na(deb_seq)] <- max(!is.na(deb_seq))
	profit_seq[(T-t0)] <- fleet_ssval_T(X[(T-t0)],sat_seq[(T-t0)],T,p,F)
	losses <- L(sat_seq,deb_seq)
	values <- as.data.frame(cbind(times,X,sat_seq,deb_seq,profit_seq,discounted_profit_seq,fleet_npv_path,losses,p,F))
	colnames(values) <- c("time","launches","satellites","debris","fleet_flowv","fleet_pv","fleet_vfn_path","collision_rate","returns","costs")
	return(values)
}
