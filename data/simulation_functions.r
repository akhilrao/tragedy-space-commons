##### functions to be used in deterministic satellite-debris model dynamic programming simulation script

# Span norm
span_norm <- function(input_matrix) {
	max(input_matrix) - min(input_matrix)
}

# Plot policy and value functions on a given grid
plot_pfn_vfn <- function(vfn,launch_pfn,removal_pfn,basegrid,labels) {
	fv_mat <- t(matrix(vfn,nrow=length(basegrid)))
	l_po_mat <- t(matrix(launch_pfn,nrow=length(basegrid)))
	r_po_mat <- t(matrix(removal_pfn,nrow=length(basegrid)))

	image2D(z=fv_mat,x=basegrid,y=basegrid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=labels[1])
	image2D(z=l_po_mat,x=basegrid,y=basegrid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100), main=labels[2])
	image2D(z=r_po_mat,x=basegrid,y=basegrid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100), main=labels[3])
}

# Build a square grid with specified curvature
build_grid <- function(gridmin, gridmax, gridlength, curvature) {
	base_piece <- (seq(from=gridmin, to=gridmax^(1/curvature), length.out=gridlength))^curvature
	sats <- rep(base_piece,times=gridlength)
	debs <- sort(rep(base_piece,times=gridlength))
	igrid <- as.data.frame(cbind(sats,debs))
	return(list(base_piece=base_piece,igrid=igrid))
}

# profit function with value function lookup
profit <- function(X,igrid,entry_no,value_fn,...) {
	S <- igrid$sats[entry_no]
	D <- igrid$debs[entry_no]
	S_next <- S_(X,S,D)
	D_next <- D_(X,S,D)
	next_state <- c(S_next,D_next)
	#print(paste0("State is (", X, "," ,S, "," ,D,")"))
	interpolation <- interpolate(next_state,igrid,value_fn)
	prof <- one_p_return(X,S) + discount_fac*interpolation
	if(is.infinite(prof)) {prof <- 0}
	return(prof)
}

# profit function with value function lookup and removal (open access)
profit_rem <- function(X,R,igrid,entry_no,value_fn,...) {
	S <- igrid$sats[entry_no]
	D <- igrid$debs[entry_no] - R
	S_next <- S_(X,S,D)
	D_next <- D_(X,S,D)
	next_state <- c(S_next,D_next)
	#print(paste0("State is (", X, "," ,S, "," ,D,")"))
	interpolation <- interpolate(next_state,igrid,value_fn)
	prof <- one_p_return_rem(X,R,S) + discount_fac*interpolation
	if(is.infinite(prof)) {prof <- 0}
	return(prof)
}

# profit function with value function lookup and removal (planner). Since it doesn't seem like I can specify a separate box constraint for removal in the optim() call, I impose it here.
profit_rem_plan <- function(policy,igrid,entry_no,value_fn,...) {
	R <- policy[2]
	if(R>igrid$debs[entry_no]) {R <- igrid$debs[entry_no]}
	S <- igrid$sats[entry_no]
	D <- igrid$debs[entry_no] - R
	S_next <- S_(policy[1],S,D)
	D_next <- D_(policy[1],S,D)
	next_state <- c(S_next,D_next)
	interpolation <- interpolate(next_state,igrid,value_fn)
	prof <- one_p_return_rem(policy[1],policy[2],S) + discount_fac*interpolation # i use the originally-chosen removal rate here so that the profit function reflects the cost of going higher, with no benefit because R is set above. This should induce the optimizer to not set policy[2] higher than the amount of debris, allowing me to pull the constraint-respecting removal rate from the optimizer's output.
	if(is.infinite(prof)) {prof <- 0}
	if(is.na(prof)) {prof <- 0}
	return(prof)
}

# Fleet's finite horizon continuation value under the given policy function - Matrix version
mat_W_ih <- function(X,V_guess,igrid,T,...) {
	n_grid_points <- sqrt(length(X))
	base_piece <- unique(igrid$sats)
	Svals <- matrix(igrid$sats,nrow=n_grid_points,ncol=n_grid_points)
	Dvals <- matrix(igrid$debs,nrow=n_grid_points,ncol=n_grid_points)
	Xvals <- matrix(X,nrow=n_grid_points,ncol=n_grid_points)
	contval <- matrix(V_guess,nrow=n_grid_points,ncol=n_grid_points)
	S_next <- S_(Xvals,Svals,Dvals)
	D_next <- D_(Xvals,Svals,Dvals)
	next_state <- c(S_next,D_next)

	for(i in 2:T) {
		interpolation <- akima::interp(x=Svals,y=Dvals,z=contval,xo=S_next,yo=D_next,linear=FALSE,extrap=TRUE)[[3]]
		contval <- one_p_return(Xvals,Svals) + discount_fac*interpolation
	}
	return(contval)
}

# Fleet's finite horizon continuation value under the given policy function
W_ih <- function(X,igrid,T,...) {
	Svals <- matrix(0,nrow=T,ncol=length(igrid$sats))
	Dvals <- matrix(0,nrow=T,ncol=length(igrid$debs))
	Xvals <- matrix(0,nrow=T,ncol=length(X))
	contval <- matrix(0,nrow=T,ncol=length(igrid$sats))

	sgl <- sqrt(dim(igrid)[1])				# "single_grid_length": assumes supplied grids are square and are stacked row-wise as columns

	Svals[1,] <- igrid$sats
	Dvals[1,] <- igrid$debs
	Xvals[1,] <- X
	contval[1,] <- one_p_return(Xvals[1,],Svals[1,])

	contval <- foreach(i=2:T, .export=ls(), .combine=rbind, .inorder=TRUE) %do% {
		Svals[i,] <- S_(Xvals[i-1,],Svals[i-1,],Dvals[i-1,])
		Dvals[i,] <- D_(Xvals[i-1,],Svals[i-1,],Dvals[i-1,])
		next_state_raw <- c(Svals[i,],Dvals[i,])
		itm <- proc.time()[3]
		next_state <- cbind(next_state_raw[1:length(X)],next_state_raw[(length(X)+1):length(next_state_raw)])
		Xvals[i,] <- foreach(j=1:length(X), .export=ls(), .combine=cbind, .inorder=TRUE) %dopar% {
			interpolate(c(next_state[j],next_state[j+length(X)]),igrid,X)
		}
		il.time <- proc.time()[3] - itm
		#print(paste0("Finished period ",i, " inner loop. Time taken: ", round(il.time,3)," seconds."))
		one_p_return(Xvals[i,],Svals[i,])*discount_fac^(i-1)
	}
	if(dim(contval)[1] < dim(Xvals)[1])  {contval <- rbind(one_p_return(Xvals[1,],Svals[1,]),contval)}
	#contval[T,] <- ((discount_fac^T)/(1-discount_fac*(1-L(Svals[T,],Dvals[T,]))))*one_p_return(L(Svals[T,],Dvals[T,])*Svals[T,],Svals[T,])
	contval[T,] <- fleet_ssval_T(Svals[T,],Dvals[T,],T)
	cv <- colSums(contval)
	return(cv)
}

# Planner's finite horizon value given a launch path
fhvf <- function(X,S,D,T,...) {
	values <- matrix(0,nrow=(T+1),ncol=5)
	values[,1] <- seq(from=0,to=T,by=1)
	values[,2] <- X
	values[1,3] <- S
	values[1,4] <- D
	values[1,5] <- p*values[1,3] - F*values[1,2]

	for(j in 2:(T+1)) {
		values[j,3] <- S_(values[(j-1),2],values[(j-1),3],values[(j-1),4])
		values[j,4] <- D_(values[(j-1),2],values[(j-1),3],values[(j-1),4])
		values[j,5] <- one_p_return(values[j,2],values[j,3])*(discount_fac^(values[(j-1),1]))
	}
	values[which(values[,4]=="NaN"),4] <- 1000000
	#values[(T+1),5] <- one_p_return(L(values[(T+1),3],values[(T+1),4])*values[(T+1),3],values[(T+1),3])*(discount_fac^(values[(T+1),1]))/(1-discount_fac*(1-L(values[(T+1),3],values[(T+1),4]))) # make final value into "steady state" value
	values[(T+1),5] <- fleet_ssval_T(values[(T+1),3],values[(T+1),4],T+1)
	fleet_npv <- sum(values[,5])
	return(fleet_npv)
}

# Planner's finite horizon value with removal
fhvf_rem <- function(XR,S,D,T,remdate,...) {
	values <- matrix(0,nrow=(T+1),ncol=6)
	values[,1] <- seq(from=0,to=T,by=1)
	values[,2] <- XR[1:(T+1)]
	values[1,3] <- S
	values[1,4] <- D
	values[,5] <- XR[(T+2):(length(XR)-1)]
	values[1,6] <- p*values[1,3] - removal_c*values[1,5] - F*values[1,2]

	for(j in 2:(T+1)) {
		if(j<remdate) {values[(j-1),5] <- 0}
		if((values[(j-1),4]<values[(j-1),5])) {values[(j-1),5] <- values[(j-1),4]}
		values[j,3] <- S_(values[(j-1),2],values[(j-1),3],(values[(j-1),4]-values[(j-1),5]) )
		values[j,4] <- D_(values[(j-1),2],values[(j-1),3],(values[(j-1),4]-values[(j-1),5]) )
		values[j,6] <- one_p_return_rem(values[j,2],values[j,5],values[j,3])*(discount_fac^(values[(j-1),1]))
	}
	values[which(values[,4]=="NaN"),4] <- 1000000
	#values[(T+1),6] <- fleet_ssval_T_rem(values[(T+1),5],values[(T+1),3],values[(T+1),4],T+1)
	fleet_npv <- sum(values[,6])

	print(values[,5])

	return(fleet_npv)
}

# Generate a series from a given policy and state grid
seriesgen_pfn <- function(X,S,D,igrid,T,...) {
	values <- matrix(0,nrow=(T+1),ncol=7)
	values[,1] <- seq(from=0,to=T,by=1)
	values[1,3] <- S
	values[1,4] <- D
	policy <- X
	
	for(i in 2:(T+1)) {
		values[(i-1),2] <- interpolate( c(values[(i-1),3],values[(i-1),4]),igrid,policy)
		values[i,3] <- S_(values[(i-1),2],values[(i-1),3],values[(i-1),4])
		values[i,4] <- D_(values[(i-1),2],values[(i-1),3],values[(i-1),4])
		values[i,5] <- one_p_return(values[(i-1),2],values[i,3])*discount_fac^(i-1)
	}
	values[(T+1),2] <- interpolate(values[T,3],igrid,policy)
	values[(T+1),5] <- fleet_ssval_T(values[(T+1),3],values[(T+1),4],T+1)
	values[which(values[,4]=="NaN"),4] <- 1000000
	values[,6] <- V_ss(values[,3],values[,4])
	values[,7] <- L(S_(values[,2],values[,3],values[,4]),D_(values[,2],values[,3],values[,4]))
	values <- as.data.frame(values)
	colnames(values) <- c("time","launches","satellites","debris","fleet_pv","satellite_pv","collision_rate")
	return(values)
}

# Generate a series from a given policy path and initial condition
seriesgen_ts <- function(X,S,D,T,...) {
	values <- matrix(0,nrow=(T+1),ncol=7)
	values[,1] <- seq(from=0,to=T,by=1)
	values[,2] <- X
	values[1,3] <- S
	values[1,4] <- D
	values[1,5] <- one_p_return(values[1,2],values[1,3])

	for(i in 2:(T+1)) {
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

# Generate a series from a given policy path and initial condition
seriesgen_ts_rem <- function(X,S,D,T,...) {
	values <- matrix(0,nrow=(T+1),ncol=8)
	values[,1] <- seq(from=0,to=T,by=1)
	values[,2] <- X[1:(T+1)]
	values[,3] <- X[(T+2):(length(X)-1)]
	values[1,4] <- S
	values[1,5] <- D
	values[1,6] <- one_p_return_rem(values[1,2],values[1,3],values[1,4])

	for(i in 2:(T+1)) {
		values[i,4] <- S_(values[(i-1),2],values[(i-1),4], (values[(i-1),5]-values[(i-1),3]) )
		values[i,5] <- D_(values[(i-1),2],values[(i-1),4], (values[(i-1),5]-values[(i-1),3]) )
		values[i,6] <- one_p_return_rem(values[i,2],values[i,3],values[i,4])*discount_fac^(i-1)
	}
	#values[(T+1),6] <- fleet_ssval_T_rem(values[(T+1),3],values[(T+1),4],values[(T+1),5],T+1)
	values[which(values[,4]=="NaN"),4] <- 1000000
	values[,7] <- V_ss(values[,4],values[,5])
	values[,8] <- L(values[,4],values[,5])
	values <- as.data.frame(values)
	colnames(values) <- c("time","launches","removals","satellites","debris","fleet_pv","satellite_pv","collision_rate")
	return(values)
}

# 2D interpolation function: interpolate the supplied values between the two grid points nearest the target in each dimension
interpolate <- function(target,grid,values) {
	# if(target[1]<max(grid$sats)&&target[2]<max(grid$debs)) {
		# print("target is in grid")
		# calculate distances from target to grid nodes
		#sgl <- sqrt(dim(grid)[1])				# "single_grid_length": assumes supplied grids are square and are stacked row-wise as columns
		# the following lines allow for rectangular grids - more points in one variable than another
		sat_grid_length <- length(unique(grid$sats))
		deb_grid_length <- length(unique(grid$debs))

		## S grid
		S_node_distances <- abs(target[1]-grid$sats)
		#S_nearest_node_idx <- sort(S_node_distances,index.return=TRUE)$ix[seq(1,2*sgl)]
		S_nearest_node_idx <- sort(S_node_distances,index.return=TRUE)$ix[seq(1,2*deb_grid_length)]
		## D_grid
		D_node_distances <- abs(target[2]-grid$debs)
		#D_nearest_node_idx <- sort(D_node_distances,index.return=TRUE)$ix[seq(1,2*sgl)]
		D_nearest_node_idx <- sort(D_node_distances,index.return=TRUE)$ix[seq(1,2*sat_grid_length)]
		
		# locate nearest S and D
		#S_nearest_nodes <- grid$sats[c(S_nearest_node_idx[1],S_nearest_node_idx[sgl+1])]
		#nearest_S_dist <- S_node_distances[c(S_nearest_node_idx[1],S_nearest_node_idx[sgl+1])]
		S_nearest_nodes <- grid$sats[c(S_nearest_node_idx[1],S_nearest_node_idx[deb_grid_length+1])]
		nearest_S_dist <- S_node_distances[c(S_nearest_node_idx[1],S_nearest_node_idx[deb_grid_length+1])]

		#D_nearest_nodes <- grid$debs[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sgl+1]))]
		#nearest_D_dist <- D_node_distances[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sgl+1]))]
		D_nearest_nodes <- grid$debs[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sat_grid_length+1]))]
		nearest_D_dist <- D_node_distances[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sat_grid_length+1]))]
		
		# "new" scheme: seems to cause problems for off-grid targets
		# calculate areas attached to each node
		# areas <- c(0,0,0,0)
		# areas[1] <- (nearest_S_dist[1]*nearest_D_dist[1])/(diff(S_nearest_nodes)*diff(D_nearest_nodes))
		# areas[2] <- (nearest_S_dist[1]*nearest_D_dist[2])/(diff(S_nearest_nodes)*diff(D_nearest_nodes))
		# areas[3] <- (nearest_S_dist[2]*nearest_D_dist[1])/(diff(S_nearest_nodes)*diff(D_nearest_nodes))
		# areas[4] <- (nearest_S_dist[2]*nearest_D_dist[2])/(diff(S_nearest_nodes)*diff(D_nearest_nodes))
		# # sanitize/normalize areas: make sure they're nonnegative and sum to 1
		# areas <- abs(areas)
		# areas <- areas/sum(areas)
		# # interpolate.
		# nearby_values <- values[intersect(S_nearest_node_idx,D_nearest_node_idx)]
		# interpolation <- nearby_values[1]*areas[4] + nearby_values[2]*areas[3] + nearby_values[3]*areas[2] + nearby_values[4]*areas[1] 
		#print(target)
		#print(c(nearest_S_dist,S_nearest_nodes,nearest_D_dist,D_nearest_nodes,areas,sum(areas)))

		# "old" scheme: seemed to work ok...?
		S_node_wts <- c(0,0)
		D_node_wts <- c(0,0)
		S_node_wts[1] <- nearest_S_dist[2]/(nearest_S_dist[1]+nearest_S_dist[2])
		S_node_wts[2] <- nearest_S_dist[1]/(nearest_S_dist[1]+nearest_S_dist[2])
		D_node_wts[1] <- nearest_D_dist[2]/(nearest_D_dist[1]+nearest_D_dist[2])
		D_node_wts[2] <- nearest_D_dist[1]/(nearest_D_dist[1]+nearest_D_dist[2])
		nearby_values <- values[intersect(S_nearest_node_idx,D_nearest_node_idx)]
		S_comp_1 <- S_node_wts[1]*nearby_values[1] + S_node_wts[2]*nearby_values[3]
		S_comp_2 <- S_node_wts[1]*nearby_values[2] + S_node_wts[2]*nearby_values[4]
		interpolation <- D_node_wts[1]*S_comp_1 + D_node_wts[2]*S_comp_2

		return(interpolation)
	# }
	#  else if(target[1]>max(grid$sats)||target[2]>max(grid$debs)) {
	#  	# print("target is off-grid")
	#  	# print(target)
	#  	return(0)
	# }
}

stock_tax_deviation <- function(tax_path,fp_launches,S,D,T,fe_eqm,...) {
	oa_launches <- oa_tsgen_stock(S,D,T,fe_eqm,tax_path)$launches[1:((T-25))]
	deviation <- sum( (oa_launches-fp_launches[1:(T-25)])^2 )
	print(deviation)
	ggplot(data=data.frame( time=seq(from=1,by=1,length=((T-25))),tax_path=tax_path )) + geom_line(aes(x=time,y=tax_path))
	return(deviation)
}

flow_tax_deviation <- function(tax_path,fp_launches,S,D,T,fe_eqm,...) {
	oa_launches <- oa_tsgen_flow(S,D,T,fe_eqm,tax_path)$launches[1:(T-25)]
	deviation <- sum( (oa_launches-fp_launches[1:(T-25)])^2 )
	print(deviation)
	ggplot(data=data.frame( time=seq(from=1,by=1,length=((T-25))),tax_path=tax_path )) + geom_line(aes(x=time,y=tax_path))
	return(deviation)
}