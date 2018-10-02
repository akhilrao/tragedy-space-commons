##### algorithms to be used in deterministic satellite-debris model dynamic programming simulation script

library(rootSolve)

### open access policy algorithm - full serial (still fast - just rootfinding with no integration)
oapolicy <- function(igrid,...) {
ptm <- proc.time()
	removals <- rep(0,length=dim(igrid)[1])
	launches <- rep(0,length=dim(igrid)[1])

	# for(j in 1:dim(igrid)[1]) {
	# 	D_temp <- igrid$debs[j]
	# 	S_temp <- igrid$sats[j]
	# 	full_rem <- ifelse( S_temp>0&&D_temp>0, D_temp/S_temp, 0.01)
	# 	Ri <- optim(par=full_rem/2, fn=satval_rem, D=D_temp, S=S_temp , launch_cost=F, removal_cost=removal_cost, method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))$par
	# 	removals[j] <- ifelse(length(Ri)==0,0,Ri*igrid$sats[j])
	# }

	for(j in 1:dim(igrid)[1]) {
		X <- uniroot.all(eqmcond,c(0,1000),S=igrid$sats[j],D=(igrid$debs[j]-removals[j]),fe_eqm=fe_eqm)
		launches[j] <- ifelse(length(X)==0,0,X)
	}

	launches[which(launches=="Inf")] <- -1
	fleet_size <- S_(launches,igrid$sats,igrid$debs)
	ih_sat_value <- V_ss(S_(launches,igrid$sats,igrid$debs),D_(launches,igrid$sats,igrid$debs))
	
	#registerDoParallel(cores=ncores)
	#ih_fleet_value <- W_ih(X=launches,igrid=igrid,T=prop_limit)
	#stopImplicitCluster()
	ih_fleet_value <- W_ss(S_(launches,igrid$sats,igrid$debs),D_(launches,igrid$sats,igrid$debs))

	loss <- L(S_(launches,igrid$sats,igrid$debs),D_(launches,igrid$sats,igrid$debs))
	satellite_oc <- rep(F*(1+r),length=length(launches))
	results <- as.data.frame(cbind(igrid,launches,removals,fleet_size,loss,ih_sat_value,satellite_oc,ih_fleet_value))
	colnames(results) <- c("satellites","debris","launches","removals","fleet_size","loss_next","satellite_value","satellite_opportunity_cost","fleet_value")
	#write.csv(results,file="./oa.policy.csv")
	print(paste("Total time taken for open access policy: ",round(((proc.time() - ptm)[3])/60,4)," minutes",sep=""))
	return(results)
}

### Computes backwards induction from given terminal value with given policy to get value function
oavalue <- function(igrid,launch_policy,removal_policy,value_fn,T) {
	count <- 0
	tot.time <- proc.time()[3]
	## make progress bar
	BIpb <- progress_bar$new(format="Doing backwards induction [:bar] :percent",total=(T+1))
	BIpb$tick(0)
	while(count<=T) {
		value_fn <- foreach(i=1:length(launch_policy), .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			result <- profit_rem(X=launch_policy[i],R=removal_policy[i],igrid=igrid,entry_no=i,value_fn=value_fn)
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

### fleet planner VFI algorithm
vfi_solver <- function(vguess,launch_pguess,removal_pguess,schedule,ncores,curvature,...) {
distlist <- list()
dev.new()
dev.new(width=10,height=10, unit="in")
par(mfrow=c(4,3))
for(g in 1:nrow(schedule)) {
	# build grid
	gridlist <- build_grid(gridmin, gridmax, schedule[g,2], curvature)
	base_piece <- gridlist[[1]]
	igrid <- gridlist[[2]]

	if(g>1&&(schedule[g,2]!=schedule[(g-1),2]) ) {
		launch_pguess <- akima::interp(x=panel$S,y=panel$D,z=panel$X,xo=base_piece,yo=base_piece,linear=TRUE,extrap=FALSE)[[3]]
		removal_pguess <- akima::interp(x=panel$S,y=panel$D,z=panel$R,xo=base_piece,yo=base_piece,linear=TRUE,extrap=FALSE)[[3]]
		vguess <- akima::interp(x=panel$S,y=panel$D,z=panel$V,xo=base_piece,yo=base_piece,linear=TRUE,extrap=FALSE)[[3]]
	}

	# initialize empty matrix for panel
	pancols <- 5
	panrows <- dim(igrid)[1]
	panel <- data.frame(matrix(0,nrow=panrows,ncol=pancols))
	colnames(panel)[1:5] <- c("V","S","D","X","R")

	# satellite-debris grid parameters and assignments for panel
	basegrid <- unique(igrid[,1])
	n_grid_points <- schedule[g,2]
	gridspace <- abs(mean(diff(igrid$sats)))
	panel$S <- igrid$sats
	panel$D <- igrid$debs

	panel$X <- as.vector(launch_pguess)
	panel$R <- as.vector(removal_pguess)
	panel$V <- as.vector(vguess)

	newX <- rep(-1,length=panrows)
	result <- cbind(newX,newX,newX,newX)

	policy.iteration.steps <- schedule[g,1]

	# initialize distance-time matrix, epsilon-delta, and count
	distmat <- data.frame(matrix(-1,nrow=10000,ncol=6))
	colnames(distmat) <- c("iteration","distance","maximization.time","policyiter.time","write.time","total.time")
	if(g<nrow(schedule)) {epsilon <- 1e-5	} #(g^0.15)*gridspace # gridspace pulls epsilon down on bigger grids, sqrt factor slows that. allow sloppiness for guess-grids
	if(g==nrow(schedule)) {epsilon <- 1e-5} # use a fixed epsilon for the final grid
	delta <- epsilon + 10
	delta2 <- delta
	count <- 0
	cat(paste("\n Stopping criteria: distance < ", epsilon, "\n", sep=""))

	# Plot initial guesses
	plot_pfn_vfn(vguess,launch_pguess,removal_pguess,basegrid,labels=c("Value fn guess","Launch policy guess","Removal policy guess"))

	# create cluster to parallelize maximization step
	registerDoParallel(cores=ncores) # do it this way instead of through makeCluster so that ls() gets passed properly

	# solver loop
	while(delta > epsilon && delta2 > epsilon) {
		## make progress bar
		#pb <- progress_bar$new(format="\nMaximizing... [:bar] :percent. Elapsed time: :elapsedfull \n",total=panrows)
		#pb$tick(0)

		## plot current value and policy functions
		plot_pfn_vfn(panel$V,panel$X,panel$R,basegrid,labels=c("Initial value fn","Initial launch policy","Initial removal policy"))

		t.tm <- proc.time()
		## maximization step
		cat(paste0("\n Starting maximization step, grid size = ", schedule[g,2]^2,": \n",". r_c: ",r_c,"\n "))
		m.tm <- proc.time()
		result <- foreach(k=1:panrows, .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			ctm <- proc.time()[3]
			solution <- optim(par=c(panel$X[k],panel$R[k]), fn=profit_rem_plan, igrid=igrid, entry_no=k, value_fn=panel$V, control=list(fnscale=-1), method="L-BFGS-B", lower=0, upper=gridmax )
			launch_rate <- solution$par[1]
			if(solution$par[2]>igrid$debs[k]) {removal_rate <- igrid$debs[k]}
			if(solution$par[2]<=igrid$debs[k]) {removal_rate <- solution$par[2]}
			vfn <- solution$value
			clock <- proc.time()[3] - ctm
			#pb$tick()
			result <- c(launch_rate,removal_rate,vfn,clock)
			return(result)
		}
		maxim.time <- (proc.time()-m.tm)[3]
		newX <- result[,1]
		newR <- result[,2]
		avg_cell_time <- mean(result[,4])
		cat(paste0("Finished. \n Average cell time: ",round(avg_cell_time,3)," seconds.\nTotal grid time: ",round((avg_cell_time*schedule[g,2]^2),3),"seconds \n"))
		rownames(newX) <- NULL

		## terminate progress bar
		#pb$terminate()

		## update panel with new policy, save old policy
		oldX <- panel$X
		oldR <- panel$R		
		panel$X <- newX
		panel$R <- newR
		newV <- result[,3]

		#cat(paste0("\n Starting propagation step... "))
		avg_period_time <- proc.time()[3]
		## calculate new value function and distance - expectation/propagation step and |V-newV| (or |XR-newXR|) 
		if(policy.iteration.steps==0 || delta>=1 ) {
			policyiter.time <- 0
			delta <- max(abs((panel$V-newV))) #infinity norm of distance between the value vectors
		}
		if(policy.iteration.steps!=0 && delta<1) {
			pit.tm <- proc.time()[3]
			policy.iteration.steps <- policy.iteration.steps
			newV <- oavalue(igrid,panel$X,panel$R,panel$V,T=policy.iteration.steps)
			policyiter.time <- proc.time()[3] - pit.tm
			delta <- max(abs((oldX-newX))) + max(abs((oldR-newR))) #sum of infinity norms of distances between the policy vectors
		}


		avg_period_time <- ifelse(policy.iteration.steps>0, (proc.time()[3] - avg_period_time)/policy.iteration.steps, 0)
		#cat(paste0("Finished. \n Average period time: ",round(avg_period_time,3)," seconds\n"))

		## update count
		iteration <- count
		count <- count+1
		## calculate |delta[i] - delta[i-1]|
		if(count>1) {delta2 <- abs(distmat$distance[count]-distmat$distance[(count-1)])}
		## update value function
		panel$V <- newV

		## plot new value function
		plot_pfn_vfn(panel$V,panel$X,panel$R,basegrid,labels=c("New value fn","New launch policy","New removal policy"))
		dev.set(dev.prev())

		## write out to csv
		wr.tm <- proc.time()
		valuemat <- matrix(panel$V,nrow=sqrt(dim(igrid)[1]),ncol=sqrt(dim(igrid)[1]),byrow=TRUE)
		launch_policymat <- matrix(panel$X,nrow=sqrt(dim(igrid)[1]),ncol=sqrt(dim(igrid)[1]),byrow=TRUE)
		removal_policymat <- matrix(panel$R,nrow=sqrt(dim(igrid)[1]),ncol=sqrt(dim(igrid)[1]),byrow=TRUE)

		vfi_vfn_filename <- paste0("nor_vfi_vfn.csv")
		vfi_lp_filename <- paste0("nor_vfi_launch_pfn.csv")
		vfi_rp_filename <- paste0("nor_vfi_removal_pfn.csv")

		write.csv(valuemat,file=vfi_vfn_filename)
		write.csv(launch_policymat,file=vfi_lp_filename)
		write.csv(removal_policymat,file=vfi_rp_filename)
		write.time <- (proc.time() - wr.tm)[3]
		## calculate total time
		tot.time <- (proc.time() - t.tm)[3]
		## fill distance-time matrix
		distmat[count,] <- c(iteration,delta,maxim.time,policyiter.time,write.time,tot.time)
		distmat_entered <- distmat[which(distmat$iteration!=-1),]
		distmatplot <- ggplot(data=distmat_entered) + 
						geom_line(aes(x=iteration,y=log(distance)),size=1) + 
						geom_hline(yintercept=log(epsilon), linetype="dashed",size=0.9) + 
						theme_minimal()
		print(distmatplot)
		dev.set(dev.next())
		## out
		cat(paste("\n Finished iteration ", count,", distance is ", delta, sep=""))
	}

	cat(paste(
		"\n Total maximization time: ", sum(distmat$maximization.time[which(distmat$maximization.time>0)]),
		"\n Total policy iteration time: ", sum(distmat$policyiter.time[which(distmat$policyiter.time>0)]),
		",\n Total CSV write time: ", sum(distmat$write.time[which(distmat$write.time>0)]),
		",\n Total algorithm time: ", sum(distmat$total.time[which(distmat$total.time>0)]), 
		",\n Final distance: ", delta, 
		",\n Total number of iterations: ", count, 
		",\n Average time per iteration: ", (sum(distmat$total.time[which(distmat$total.time>0)]))/count,"\n\n", sep=""))
	distlist[[g]] <- cbind(distmat[which(distmat$iteration!=-1),],g,schedule[g,2])
}
	stopImplicitCluster()

	distlistvec <- as.data.frame(rbindlist(distlist))
	distlistvecplot <- ggplot(data=distlistvec) + 
						geom_line(aes(x=iteration,y=log(distance), group=g, color=factor(g)),size=1) + 
						scale_color_viridis(option="inferno",discrete=TRUE) +
						theme_minimal()
	#print(distlistvecplot)

	return(as.data.frame(cbind(satellites=panel$S,debris=panel$D,optimal_fleet_vfn=panel$V,optimal_launch_pfn=panel$X,optimal_removal_pfn=panel$R,optimal_fleet_size=S_(panel$X,panel$S,panel$D),optimal_sat_val=V_ss(S_(panel$X,panel$S,panel$D),panel$S,panel$D))))
	dev.off()
	dev.off()
}

# fleet planner time-series value function approximator
fpts_vfsolve <- function(igrid,T,ncores,...) {
	registerDoParallel(cores=ncores)

	policy_value <- as.data.frame(matrix(-1,nrow=dim(igrid)[1],ncol=6))

	print("Generating finite-horizon guess...")
	policy_value <- foreach(j=1:dim(igrid)[1], .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
		ctm <- proc.time()[3]
		parms <- c(rep(0.01,length=T),0,igrid$sats[j],igrid$debs[j])
		solution <- optim(par=parms[1:(T+1)],fn=fhvf, S=parms[(T+2)], D=parms[(T+3)], T=T, control=list(fnscale=-1, pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)
		policy <- solution$par[1]
		value <- solution$value
		satval <- V_ss(S_(policy,igrid$sats,igrid$debs),igrid$debs)[1]
		clock <- proc.time()[3] - ctm
		result <- c(igrid$sats[j],igrid$debs[j],policy,satval,value,clock)
		result
	} 
	policy_value <- as.data.frame(policy_value)
	colnames(policy_value) <- c("satellites","debris","launches","sat_value","fleet_value","compute_time")
	avg_time_per_cell <- mean(policy_value$compute_time)
	total_time <- sum(policy_value$compute_time)/(60*60*32)
	print(paste0("Finished. Average time per cell: ", round(avg_time_per_cell,3)," seconds. Total wall time: ", round(total_time,3)," hours."))

	stopImplicitCluster()

	return(policy_value)
}

# fleet planner time series generator
fp_tsgen <- function(S,D,T,...) {
	ctm <- proc.time()[3]
	parms <- c(rep(0.01,length=T),0,S,D)
	launch_path <- optim(par=parms[1:(T+1)],fn=fhvf, S=parms[(T+2)], D=parms[(T+3)], T=T, control=list(fnscale=-1, pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)$par
	launch_path_clock <- proc.time()[3] - ctm

	ctm <- proc.time()[3]
	time_series <- seriesgen_ts(launch_path,S,D,T)
	propagation_clock <- proc.time()[3] - ctm

	print("Fleet planner time series generated.") 
	print(paste0("Time to compute launch path: ", launch_path_clock, " seconds."))
	print(paste0("Time to propagate stocks: ", propagation_clock, " seconds."))

	return(time_series)
}

# Generate an optimal time series with endogenous removal
fp_tsgen_endorem<- function(S,D,T,fe_eqm,removal_c,remstart,...) {
	ctm <- proc.time()[3]
	parms <- c(rep(0.01,length= (2*T+3) ),S,D)
	launch_path <- optim(par=parms[1:(2*T+3)],fn=fhvf_rem, S=parms[(2*T+4)], D=parms[(2*T+5)], T=T, remdate=remstart, control=list(fnscale=-1, pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)$par
	launch_path_clock <- proc.time()[3] - ctm

	ctm <- proc.time()[3]
	time_series <- seriesgen_ts_rem(launch_path,S,D,T)
	propagation_clock <- proc.time()[3] - ctm

	print("Fleet planner time series generated.") 
	print(paste0("Time to compute launch path: ", launch_path_clock, " seconds."))
	print(paste0("Time to propagate stocks: ", propagation_clock, " seconds."))

	return(time_series)
}

# open access time series generator
oa_tsgen <- function(S,D,T,fe_eqm,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:T ) {
		X <- uniroot.all(eqmcond,c(0,100000),S=sat_path[t],D=D,fe_eqm=fe_eqm)
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
		fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

# Generate an open access time series with endogenous removal
oa_tsgen_endorem<- function(S,D,T,fe_eqm,removal_c,remstart,epsilon,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- c(S,rep(0,length=(T)))
	deb_path <- c(D,rep(0,length=(T)))
	fleet_pv_path <- rep(0,length=(T+1))
	removal_R <- rep(0,length=(T+1))
	count <- 1

	delta <- epsilon+10
	print(paste0("Beginning fixed-point iteration. Convergence criteria: delta < ", epsilon,"."))

	launch_policy <- ggplot(data=data.frame(time=seq(0,T),launches=launch_path)) + 
							geom_line(aes(x=time,y=launches))
	removal_policy <- ggplot(data=data.frame(time=seq(0,T),removals=removal_R)) + 
							geom_line(aes(x=time,y=removals))

	grid.arrange(launch_policy,removal_policy,ncol=2)
	if(remstart>1){
		for(t in 1:(remstart-1) ) {
			# before the technology: no one can use it, OA as usual
			X <- uniroot.all(eqmcond,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm)
			if(length(X)==0) { launch_path[t] <- 0 }	
			if(length(X)>0) { launch_path[t] <- X }
			sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
			deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
			fleet_pv_path[t] <- one_p_return_rem(launch_path[t],removal_R[t],sat_path[t])*discount_fac^(t-1)
		}
	}
	while(delta>=epsilon){
		print(paste0("Starting iteration ", count,"."))
		oldX <- launch_path
		oldR <- removal_R
		for(t in 1:(T-1) ) {
			
			# the day before the technology, launchers think about how much they'll buy
			if(t==remstart){
				X <- uniroot.all(eqmcond_exorem,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm,r_cbar=removal_c[t+1] ,Rbar=removal_R[t+1])
				if(length(X)==0) { launch_path[t] <- oldX[t] }	
				if(length(X)>0) { launch_path[t] <- X }
				sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
				deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
				fleet_pv_path[t] <- one_p_return_rem(launch_path[t],removal_R[t],sat_path[t])*discount_fac^(t-1)
			}
			# when the technology is out there: owners remove, launchers anticipate their ownership
			if(t>=(remstart+1) ){
				full_rem <- ifelse(sat_path[t]>0,deb_path[t]/sat_path[t],0)
				solution <- optim(par=full_rem/2, fn=satval_rem, D=deb_path[t], S=sat_path[t] , launch_cost=F, removal_cost= (removal_c[t]*F) , method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))
				removal_R[t] <- solution$par
				deb_path[t] <- deb_path[t]-(sat_path[t]*removal_R[t])
				#print(deb_path[t])
				#print(solution$value)
				X <- uniroot.all(eqmcond_exorem,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm,r_cbar=removal_c[t+1],Rbar=removal_R[t+1])
				if(length(X)==0) { launch_path[t] <- oldX[t] }	
				if(length(X)>0) { launch_path[t] <- X }
				sat_path[t+1] <- S_(launch_path[t],sat_path[t], deb_path[t] )
				deb_path[t+1] <- D_(launch_path[t],sat_path[t], deb_path[t] )
				fleet_pv_path[t] <- one_p_return_rem(launch_path[t],removal_R[t],sat_path[t])*discount_fac^(t-1)
			}
		}
		removal_R[T] <- optim(par=deb_path[T]/2, fn=satval_rem, D=deb_path[T], S=sat_path[T] , launch_cost=F, removal_cost=removal_c[T], method="L-BFGS-B", lower=0, upper=deb_path[T], control=list(fnscale=-1))$par
		#Sys.sleep(0)
		launch_policy <- launch_policy +
							geom_line(data=data.frame(time=seq(0,T),launches=launch_path), aes(x=time,y=launch_path), linetype="dashed")
		removal_policy <- removal_policy +
							geom_line(data=data.frame(time=seq(0,T),removals=removal_R), aes(x=time,y=removals), linetype="dashed")
		grid.arrange(launch_policy,removal_policy,ncol=2)

		delta <- 0.5*max(abs(oldX - launch_path)) + 0.5*max(abs(oldR - removal_R))
		print(paste0("Finished iteration ", count, ". Mean supnorm of launch and removal policies: ",delta,"."))
		count <- count + 1
	}

	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path),agg_removal=removal_R)

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}


# Generate an open access time series with exogenous removal
oa_tsgen_exorem<- function(S,D,T,fe_eqm,removal_c,removal_R,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:(T-1) ) {
		X <- uniroot.all(eqmcond_exorem,c(0,100000),S=sat_path[t],D=(deb_path[t]-removal_R[t]),fe_eqm=fe_eqm,r_cbar=removal_c[t+1],Rbar=removal_R[t+1])
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t]-removal_R[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t]-removal_R[t])
		fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

# Generate a series from a given policy model and initial condition
seriesgen_model <- function(S,D,T,coefs,igrid,degapprox,...) {
	values <- matrix(0,nrow=(T+1),ncol=7)
	values[,1] <- seq(from=0,to=T,by=1)
	values[,2] <- 0
	values[1,3] <- S
	values[1,4] <- D
	values[1,5] <- one_p_return(values[1,2],values[1,3])

	for(i in 2:(T+1)) {
		translated_state <- t(as.matrix(translate(c(values[(i-1),3],values[(i-1),4]),igrid)))
		colnames(translated_state) <- c("S","D")
		translated_state <- as.data.frame(translated_state)
		translated_state_basis <- as.matrix(basis(translated_state,degapprox,isV=0))
		X <- translated_state_basis%*%coefs
		X <- untranslate(X,igrid)
		values[(i-1),2] <- X
		values[i,3] <- S_(values[(i-1),2],values[(i-1),3],values[(i-1),4])
		values[i,4] <- D_(values[(i-1),2],values[(i-1),3],values[(i-1),4])
		values[i,5] <- one_p_return(values[(i-1),2],values[i,3])
	}
	values[(T+1),5] <- fleet_ssval_T(values[(T+1),3],values[(T+1),4],T+1)
	values[which(values[,4]=="NaN"),4] <- 1000000
	values[,6] <- V_ss(values[,3],values[,4])
	values[,7] <- L(values[,3],values[,4])
	values <- as.data.frame(values)
	colnames(values) <- c("time","launches","satellites","debris","fleet_pv","satellite_pv","collision_rate")
	return(values)
}

# Generate an open access time series with a stock control
oa_tsgen_stock <- function(S,D,T,fe_eqm,price_path,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:T ) {
		X <- uniroot.all(eqmcond_stock,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm,stock=price_path[t+1])
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
		fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

# Generate an open access time series with a flow control
oa_tsgen_flow <- function(S,D,T,fe_eqm,price_path,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:(T-1) ) {
		X <- uniroot.all(eqmcond_flow,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm,flow_t=price_path[t],flow_t1=price_path[t+1])
		if(length(X)==0) { launch_path[t] <- 0 }	
		if(length(X)>0) { launch_path[t] <- X }
		sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
		deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
		fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

oa_tsgen_flow_shutdown <- function(S,D,T,fe_eqm,price_path,shutdown_date,...) {
	ctm <- proc.time()[3]
	launch_path <- rep(0,length=(T+1))
	sat_path <- rep(S,length=(T+1))
	deb_path <- rep(D,length=(T+1))
	fleet_pv_path <- rep(0,length=(T+1))
	for(t in 1:(T-1) ) {
		if(t<shutdown_date) {
			X <- uniroot.all(eqmcond_flow,c(0,100000),S=sat_path[t],D=deb_path[t],fe_eqm=fe_eqm,flow_t=price_path[t],flow_t1=price_path[t+1])
			if(length(X)==0) { launch_path[t] <- 0 }	
			if(length(X)>0) { launch_path[t] <- X }
			sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
			deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
		}
		if(t>=shutdown_date) {
			launch_path[t] <- 0
			sat_path[t+1] <- S_(launch_path[t],sat_path[t],deb_path[t])
			deb_path[t+1] <- D_(launch_path[t],sat_path[t],deb_path[t])
		}
		fleet_pv_path[t] <- one_p_return(launch_path[t],sat_path[t])*discount_fac^(t-1)
	}
	fleet_pv_path[T+1] <- fleet_ssval_T(sat_path[T],deb_path[T],T)
	clock <- proc.time()[3] - ctm

	time_series <- data.frame(time=seq(0,T),launches=launch_path,satellites=sat_path,debris=deb_path,fleet_pv=fleet_pv_path,satellite_pv=V_ss(launch_path,sat_path),collision_rate=L(sat_path,deb_path))

	print("Open access time series generated.") 
	print(paste0("Time to compute launch path and propagate stocks: ", clock, " seconds."))

	return(time_series)
}

# Calculate optimal stock and flow control taxes by minimizing the deviation of the open access path from a pre-computed optimal path
optimal_stock_tax_path <- function(fp_launches,S,D,fe_eqm,...) {
	T <- length(fp_launches)
	tax_init <- fp_launches
	tax_path <- optim(tax_init,stock_tax_deviation,S=S,D=D,fp_launches=fp_launches, T=T, fe_eqm=fe_eqm, control=list(pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)$par
	return(tax_path)
}

optimal_flow_tax_path <- function(fp_launches,S,D,fe_eqm,...) {
	T <- length(fp_launches)
	tax_init <- fp_launches
	tax_path <- optim(tax_init,flow_tax_deviation,S=S,D=D,fp_launches=fp_launches, T=T, fe_eqm=fe_eqm, control=list(pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)$par
	return(tax_path)
}
