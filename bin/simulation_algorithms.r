##### algorithms to be used in deterministic satellite-debris model dynamic programming simulation script

library(rootSolve)

# fleet planner time series generator
fp_tsgen <- function(S,D,T,fe_eqm,launch_con,asats_seq,...) {
	ctm <- proc.time()[3]
	parms <- c(rep(15,length=T),S,D)
	launch_path <- optim(par=parms[1:T],fn=fhvf, S=parms[(T+1)], D=parms[(T+2)], asats_seq=asats_seq, launch_con=launch_con, T=T, control=list(fnscale=-1, pgtol=1e-2),method="L-BFGS-B",lower=0,upper=1000)$par
	launch_path_clock <- proc.time()[3] - ctm

	ctm <- proc.time()[3]
	time_series <- seriesgen_ts(launch_path,S,D,T,asats_seq)
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

### fleet planner VFI algorithm: given a contval, solve for a policy
dynamic_vfi_solver <- function(vguess,launch_pguess,gridmin,gridmax,t,contval,...) {
distlist <- list()
dev.new()
dev.new(width=20,height=15, unit="in")
par(mfrow=c(4,3))
	# build grid
	gridlist <- build_grid(gridmin, gridmax, 25, 1)
	base_piece <- gridlist[[1]]
	igrid <- gridlist[[2]]

	# initialize empty matrix for panel
	pancols <- 4
	panrows <- dim(igrid)[1]
	panel <- data.frame(matrix(0,nrow=panrows,ncol=pancols))
	colnames(panel)[1:4] <- c("V","S","D","X")

	# satellite-debris grid parameters and assignments for panel
	basegrid <- unique(igrid[,1])
	n_grid_points <- length(basegrid)^2
	gridspace <- abs(mean(diff(igrid$sats)))
	panel$S <- igrid$sats
	panel$D <- igrid$debs

	panel$X <- as.vector(launch_pguess)
	panel$V <- as.vector(contval)

	newX <- rep(-1,length=panrows)
	result <- cbind(newX,newX,newX)

	# initialize distance-time matrix, epsilon-delta, and count
	distmat <- data.frame(matrix(-1,nrow=1000,ncol=6))
	colnames(distmat) <- c("iteration","distance","maximization.time","policyiter.time","write.time","total.time")
	epsilon <- 1e-4
	delta <- epsilon + 10
	delta2 <- delta
	count <- 0
	cat(paste("\n Stopping criteria: distance < ", epsilon, "\n", sep=""))

	# Plot initial guesses
	plot_pfn_vfn(vguess,launch_pguess,basegrid,labels=c("Value fn guess","Launch policy guess"))

	# solver loop
	while(delta > epsilon && delta2 > epsilon) {

		## plot current value and policy functions
		plot_pfn_vfn(panel$V,panel$X,basegrid,labels=c("Initial value fn","Initial launch policy"))

		t.tm <- proc.time()
		## maximization step
		cat(paste0("\n Starting maximization step, grid size = ", n_grid_points,": \n",". r_c: ",r_c,"\n "))
		m.tm <- proc.time()
		result <- foreach(k=1:panrows, .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			ctm <- proc.time()[3]
			solution <- optim(par=panel$X[k], fn=fleet_preval, igrid=igrid, entry_no=k, value_fn=panel$V, asats=asats[t], t=t, control=list(fnscale=-1), method="L-BFGS-B", lower=gridmin, upper=gridmax)
			launch_rate <- solution$par[1]
			vfn <- solution$value
			clock <- proc.time()[3] - ctm
			result <- c(launch_rate,vfn,clock)
			return(result)
		}

		maxim.time <- (proc.time()-m.tm)[3]
		newX <- result[,1]
		avg_cell_time <- mean(result[,3])
		cat(paste0("Finished. \n Average cell time: ",round(avg_cell_time,3)," seconds.\nTotal grid time: ",round((avg_cell_time*n_grid_points),3),"seconds \n"))
		rownames(newX) <- NULL

		## update panel with new policy, save old policy
		oldX <- panel$X
		panel$X <- newX
		newV <- result[,2]

		if(delta2<1) {policy.iteration.steps <- policy.iteration.steps+25}

		#cat(paste0("\n Starting propagation step... "))
		avg_period_time <- proc.time()[3]
		## calculate new value function and distance - expectation/propagation step and |V-newV| (or |X-newX|) 
		if(policy.iteration.steps==0 || delta>=1 ) {
			policyiter.time <- 0
			delta <- max(abs((panel$V-newV)))
		}
		if(policy.iteration.steps!=0 && delta<1) {
			pit.tm <- proc.time()[3]
			policy.iteration.steps <- policy.iteration.steps
			newV <- oavalue(igrid,panel$X,panel$R,panel$V,T=policy.iteration.steps)
			policyiter.time <- proc.time()[3] - pit.tm
			delta <- max(abs((oldX-newX)))
		}

		avg_period_time <- ifelse(policy.iteration.steps>0, (proc.time()[3] - avg_period_time)/policy.iteration.steps, 0)

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

	return(as.data.frame(cbind(satellites=panel$S,debris=panel$D,optimal_fleet_vfn=panel$V,optimal_launch_pfn=panel$X,optimal_fleet_size=S_(panel$X,panel$S,panel$D))))
	dev.off()
	dev.off()
}

### fleet planner VFI algorithm: obtain steady state
ss_vfi_solver <- function(vguess,launch_pguess,gridmin,gridmax,...) {
distlist <- list()
dev.new()
dev.new(width=20,height=15, unit="in")
par(mfrow=c(4,3))
	# build grid
	gridlist <- build_grid(gridmin, gridmax, 25, 1)
	base_piece <- gridlist[[1]]
	igrid <- gridlist[[2]]

	# initialize empty matrix for panel
	pancols <- 4
	panrows <- dim(igrid)[1]
	panel <- data.frame(matrix(0,nrow=panrows,ncol=pancols))
	colnames(panel)[1:4] <- c("V","S","D","X")

	# satellite-debris grid parameters and assignments for panel
	basegrid <- unique(igrid[,1])
	n_grid_points <- length(basegrid)^2
	gridspace <- abs(mean(diff(igrid$sats)))
	panel$S <- igrid$sats
	panel$D <- igrid$debs

	panel$X <- as.vector(launch_pguess)
	panel$V <- as.vector(vguess)

	newX <- rep(-1,length=panrows)
	result <- cbind(newX,newX,newX)

	# initialize distance-time matrix, epsilon-delta, and count
	distmat <- data.frame(matrix(-1,nrow=1000,ncol=6))
	colnames(distmat) <- c("iteration","distance","maximization.time","policyiter.time","write.time","total.time")
	epsilon <- 1e-4
	delta <- epsilon + 10
	delta2 <- delta
	count <- 0
	cat(paste("\n Stopping criteria: distance < ", epsilon, "\n", sep=""))

	# Plot initial guesses
	plot_pfn_vfn(vguess,launch_pguess,basegrid,labels=c("Value fn guess","Launch policy guess"))

	# solver loop
	while(delta > epsilon && delta2 > epsilon) {

		## plot current value and policy functions
		plot_pfn_vfn(panel$V,panel$X,basegrid,labels=c("Initial value fn","Initial launch policy"))

		t.tm <- proc.time()
		## maximization step
		cat(paste0("\n Starting maximization step, grid size = ", n_grid_points,": \n",". r_c: ",r_c,"\n "))
		m.tm <- proc.time()
		result <- foreach(k=1:panrows, .export=ls(), .combine=rbind, .inorder=TRUE) %dopar% {
			ctm <- proc.time()[3]
			solution <- optim(par=panel$X[k], fn=fleet_preval, igrid=igrid, entry_no=k, value_fn=panel$V, asats=asats[length(asats)], t=length(asats), control=list(fnscale=-1), method="L-BFGS-B", lower=gridmin, upper=gridmax)
			launch_rate <- solution$par[1]
			vfn <- solution$value
			clock <- proc.time()[3] - ctm
			result <- c(launch_rate,vfn,clock)
			return(result)
		}

		maxim.time <- (proc.time()-m.tm)[3]
		newX <- result[,1]
		avg_cell_time <- mean(result[,3])
		cat(paste0("Finished. \n Average cell time: ",round(avg_cell_time,3)," seconds.\nTotal grid time: ",round((avg_cell_time*n_grid_points),3),"seconds \n"))
		rownames(newX) <- NULL

		## update panel with new policy, save old policy
		oldX <- panel$X
		panel$X <- newX
		newV <- result[,2]

		if(delta2<1) {policy.iteration.steps <- policy.iteration.steps+25}

		#cat(paste0("\n Starting propagation step... "))
		avg_period_time <- proc.time()[3]
		## calculate new value function and distance - expectation/propagation step and |V-newV| (or |X-newX|) 
		if(policy.iteration.steps==0 || delta>=1 ) {
			policyiter.time <- 0
			delta <- max(abs((panel$V-newV)))
		}
		if(policy.iteration.steps!=0 && delta<1) {
			pit.tm <- proc.time()[3]
			policy.iteration.steps <- policy.iteration.steps
			newV <- oavalue(igrid,panel$X,panel$R,panel$V,T=policy.iteration.steps)
			policyiter.time <- proc.time()[3] - pit.tm
			delta <- max(abs((oldX-newX)))
		}

		avg_period_time <- ifelse(policy.iteration.steps>0, (proc.time()[3] - avg_period_time)/policy.iteration.steps, 0)

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

	return(as.data.frame(cbind(satellites=panel$S,debris=panel$D,optimal_fleet_vfn=panel$V,optimal_launch_pfn=panel$X,optimal_fleet_size=S_(panel$X,panel$S,panel$D))))
	dev.off()
	dev.off()
}

