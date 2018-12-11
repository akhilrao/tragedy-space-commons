##### functions to be used in deterministic satellite-debris model dynamic programming simulation script

# Planner's finite horizon value given a launch path
fhvf <- function(X,S,D,T,asats_seq,launch_con,...) {
	times <- seq(from=1,to=T,by=1)	
	sat_seq <- rep(0,length=T)
	deb_seq <- rep(0,length=T)
	profit_seq <- rep(0,length=T)

	for(t in 1:length(X)) {
		if(X[t]>launch_con[t]) {X[t] <- launch_con[t]}
	}

	sat_seq[1] <- S
	deb_seq[1] <- D
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1)

	for(j in 2:T) {
		sat_seq[j] <- S_(X[(j-1)],sat_seq[(j-1)],deb_seq[(j-1)])
		deb_seq[j] <- D_(X[(j-1)],sat_seq[(j-1)],deb_seq[(j-1)],asats_seq[(j-1)])
		profit_seq[j] <- one_p_return(X[j],sat_seq[j],j)*(discount_fac^(times[(j-1)]))
		print(profit_seq[j])
	}
	deb_seq[is.na(deb_seq)] <- min(max(!is.na(deb_seq)),D)
	profit_seq[T] <- fleet_ssval_T(sat_seq[T],deb_seq[T],T)
	fleet_npv <- sum(profit_seq)
	return(fleet_npv)
}

# Generate a series from a given policy path and initial condition
seriesgen_ts <- function(X,S,D,T,asats_seq,...) {
	times <- seq(from=1,to=T,by=1)	
	sat_seq <- rep(0,length=T)
	deb_seq <- rep(0,length=T)
	profit_seq <- rep(0,length=T)
	
	sat_seq[1] <- S
	deb_seq[1] <- D
	profit_seq[1] <- one_p_return(X[1],sat_seq[1],1)

	for(k in 2:T) {
		sat_seq[k] <- S_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)])
		deb_seq[k] <- D_(X[(k-1)],sat_seq[(k-1)],deb_seq[(k-1)],asats_seq[(k-1)])
		profit_seq[k] <- one_p_return(X[k],sat_seq[k],k)*(discount_fac^(times[(k-1)]))
	}
	deb_seq[is.na(deb_seq)] <- max(!is.na(deb_seq))
	profit_seq[T] <- fleet_ssval_T(sat_seq[T],deb_seq[T],T)
	losses <- L(sat_seq,deb_seq)
	values <- as.data.frame(cbind(times,X,sat_seq,deb_seq,profit_seq,losses))
	colnames(values) <- c("time","launches","satellites","debris","fleet_pv","collision_rate")
	return(values)
}

# Plot policy and value functions on a given grid
plot_pfn_vfn <- function(vfn,launch_pfn,basegrid,labels) {
	fv_mat <- t(matrix(vfn,nrow=length(basegrid)))
	l_po_mat <- t(matrix(launch_pfn,nrow=length(basegrid)))

	image2D(z=fv_mat,x=basegrid,y=basegrid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100),contour=TRUE,main=labels[1])
	image2D(z=l_po_mat,x=basegrid,y=basegrid,xlab=c("Debris"),ylab=c("Satellites"),col=plasma(n=100), main=labels[2])
}

# Build a square grid with specified curvature
build_grid <- function(gridmin, gridmax, gridlength, curvature) {
	base_piece <- (seq(from=gridmin, to=gridmax^(1/curvature), length.out=gridlength))^curvature
	sats <- rep(base_piece,times=gridlength)
	debs <- sort(rep(base_piece,times=gridlength))
	igrid <- as.data.frame(cbind(sats,debs))
	return(list(base_piece=base_piece,igrid=igrid))
}

# Convert a square grid to a panel
grid_to_panel <- function(gridlist,launch_pguess,contval) {
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

	return(panel)
}

# 2D interpolation function: interpolate the supplied values between the two grid points nearest the target in each dimension
interpolate <- function(target,grid,values) {
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
		S_nearest_nodes <- grid$sats[c(S_nearest_node_idx[1],S_nearest_node_idx[deb_grid_length+1])]
		nearest_S_dist <- S_node_distances[c(S_nearest_node_idx[1],S_nearest_node_idx[deb_grid_length+1])]

		D_nearest_nodes <- grid$debs[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sat_grid_length+1]))]
		nearest_D_dist <- D_node_distances[sort(c(D_nearest_node_idx[1],D_nearest_node_idx[sat_grid_length+1]))]
		
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
}
