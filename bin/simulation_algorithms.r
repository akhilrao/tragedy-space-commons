##### algorithms to be used in deterministic satellite-debris model dynamic programming simulation script

library(rootSolve)

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
	# fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:T ) {
		X <- uniroot.all(eqmcond,c(0,1e+6),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm[t+1],asats=asats[t])
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		if(launch_path[t]>launch_con[t]) {launch_path[t] <- launch_con[t]}
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t],asats[t])
		# fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	# fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	# time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))
	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

### open access policy algorithm - full serial (still fast - just rootfinding with no integration)
oapolicy <- function(igrid,fe_eqm,t,...) {
ptm <- proc.time()
	launches <- rep(0,length=dim(igrid)[1])

	for(j in 1:dim(igrid)[1]) {
		X <- uniroot.all(eqmcond,c(0,1000),S=igrid$sats[j],D=(igrid$debs[j]),fe_eqm=fe_eqm[t+1],asats=asats[t])
		launches[j] <- ifelse(length(X)==0,0,X)
		if(launches[j]>launch_con[j]) {launches[j] <- launch_con[j]}
	}

	launches[which(launches=="Inf")] <- -1
	fleet_size <- S_(launches,igrid$sats,igrid$debs)

	loss <- L(S_(launches,igrid$sats,igrid$debs),D_(launches,igrid$sats,igrid$debs))
	results <- as.data.frame(cbind(igrid,launches,fleet_size,loss))
	colnames(results) <- c("satellites","debris","launches","fleet_size","loss_next")
	#write.csv(results,file="./oa.policy.csv")
	print(paste("Total time taken for open access policy: ",round(((proc.time() - ptm)[3])/60,4)," minutes",sep=""))
	return(results)
}

### Computes backwards induction from given terminal value with given policy to get value function
policy_eval_BI <- function(igrid,launch_policy,value_fn,T) {
	count <- 0
	tot.time <- proc.time()[3]
	## make progress bar
	BIpb <- progress_bar$new(format="Doing backwards induction [:bar] :percent",total=(T+1))
	BIpb$tick(0)
	while(count<=T) {
		value_fn <- foreach(i=1:length(launch_policy), .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			result <- fleet_preval(X=launch_policy[i],S=igrid$sats[i],D=igrid$debs[i],value_fn=value_fn,asats=asats,t=T,p=p,F=F,igrid=igrid)
			result
		}
		count <- count + 1
		BIpb$tick()
	}
	BIpb$terminate()
	tot.time <- round(proc.time()[3] - tot.time,3)/60
	cat(paste0("\nTotal time taken for backwards induction: ", tot.time, " minutes. Average time per period: ", tot.time/length(launch_policy), " minutes.\n"))
	return(value_fn)
}

### function to simulate fleet planner paths from different time points
opt_path_from_t <- function(pi_t,F_t,St,Dt,t,gridmin,gridmax,...) {
	T <- length(pi_t)
	policy_list <- obtain_policy_fns(gridmin,gridmax,T)
	optimal_path <- path_from_policies(pi_t,F_t,St,Dt,t)
}

### function to compute policy functions along a given returns and cost path
obtain_policy_fns <- function(gridmin,gridmax,T,...) {
	final_policy <- ss_vfi_solver(vguess,launch_pguess,gridmin,gridmax,T)
	policy_list <- list()
	policy_list[[1]] <- dynamic_vfi_solver(vguess,launch_pguess,gridmin,gridmax,1,final_policy)
	policy_list[[T]] <- final_policy
	for(i in 2:(T-1)) {
		policy_list[[i]] <- dynamic_vfi_solver(vguess,launch_pguess,gridmin,gridmax,i,policy_list[[i-1]])
	}
	return(policy_list)
}

### fleet planner VFI algorithm: solve for policy assuming steady state reached or assuming a perfect-foresight path to steady state
dynamic_vfi_solver <- function(panel,igrid,asats,t,T,p,F,...) {
	# initialize hyperparameters
	panrows <- nrow(panel)
	gridmin <- min(igrid)
	gridmax <- max(igrid)
	n_grid_points <- length(unique(igrid[,1]))^2
	base_grid <- unique(igrid[,1])
	policy.evaluation.steps <- 25
	dev.new(width=12,height=5,unit="in")
	par(mfrow=c(1,2))

	# initialize output objects
	newX <- rep(-1,length=panrows)
	result <- cbind(newX,newX,new)
	# initialize epsilon-delta and count
	ifelse(t==T, epsilon <- 1e-3, epsilon <- 2e-1) # tighter epsilon for value function convergence in final period, looser epsilon for policy function convergence in prior periods.
	ifelse(t==T,panel$X <- panel$X, panel$X <- rnorm(length(panel$X),mean=10,sd=1))
	delta_old <- 0
	delta <- epsilon + 10
	delta2 <- 25
	count <- 0
	cat(paste("\n Stopping criteria: distance < ", epsilon, "\n", sep=""))

	# solver loop
	while(delta > epsilon) {

		## plot pfn and vfn
		plot_pfn_vfn(panel$V,panel$X,base_grid,c("Value function","Policy function"))

		t.tm <- proc.time()
		## maximization step
		cat(paste0("\n Starting maximization step, grid size = ", n_grid_points,": \n","."))
		m.tm <- proc.time()
		result <- foreach(k=1:panrows, .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			ctm <- proc.time()[3]
			solution <- optim(par=panel$X[k], fn=fleet_preval, S=panel$S[k], D=panel$D[k], value_fn=panel$V, asats=asats, t=t, p=p, F=F, igrid=igrid, control=list(fnscale=-1), method="L-BFGS-B", lower=gridmin, upper=gridmax)
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
			delta <- max(abs((panel$V-newV)))
			panel$V <- newV
			panel$X <- newX
		}
		# if(t==T && delta2<0.01) {
		# 	newV <- policy_eval_BI(igrid,panel$X,panel$V,T=policy.evaluation.steps)
		# 	delta <- max(abs((panel$X-newX)))
		# 	panel$V <- newV
		# 	policy.evaluation.steps <- policy.evaluation.steps+25
		# }
		if(t!=T) {
			delta <- max(abs((panel$X-newX)))
			panel$X <- newX
		}

		## update old delta and calculate percentage change in delta by midpoint method		
		delta2 <- abs( (delta-delta_old)/(0.5*(delta+delta_old)) )
		delta_old <- delta

		## update count
		iteration <- count
		count <- count+1
	
		## calculate total time
		tot.time <- (proc.time() - t.tm)[3]
		## out
		cat(paste("\n Finished iteration ", count,", distance is ", delta, sep=""))
	}

	if(t!=T) {panel$V <- newV}

	return(as.data.frame(cbind(satellites=panel$S,debris=panel$D,optimal_launch_pfn=panel$X,optimal_fleet_vfn=panel$V,optimal_fleet_size=S_(panel$X,panel$S,panel$D),t=t,p=p[t],F=F[t])))
}
