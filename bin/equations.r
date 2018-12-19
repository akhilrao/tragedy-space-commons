##### functions to be used in deterministic satellite-debris model dynamic programming solver script

# Collision probability 
L <- function(S,D,...) {
	#1 - exp(-aSS*S-aSD*D)			#negexp rate
	#SD <- ifelse(S*D==0,1e-7,S*D)
	pmax(pmin(aS*S + aD*D + aSS*S^2 + aSD*(S*D) + aDD*D^2,1),0)		#statmech rate
	#pmax(pmin(aS*S + aD*D + aSS*S^2 + aSD*(S*D),1),0)		#statmech rate
}

# Derivatives of collision probability
L_S <- function(S,D,...) {
	ifelse(L(S,D)==1, 0, aS + 2*aSS*S + aSD*D)
}

L_D <- function(S,D,...) {
	ifelse(L(S,D)==1, 0, aD + 2*aDD*D + aSD*S)
}

# debris growth function
G <- function(S,D,...) {
	sat_caused <- ifelse( (S+D)>0, L(S,D)*(S/(S+D)), 0)
	deb_caused <- ifelse( (S+D)>0, L(S,D)*(D/(S+D)), 0)
	newfrags <- bSS*sat_caused*S + bSD*deb_caused*S + aDDbDD*D^2
	return(newfrags)
}

# Satellite law of motion 
S_ <- function(X,S,D,...) {
	S*(1-avg_sat_decay)*(1-L(S,D)) + X
}

# Debris law of motion
D_ <- function(X,S,D,asats,...) {
	stock <- D*(1-d) + G(S,D) + m*X + Z_coef*avg_sat_decay*S + asat_coef*asats
	stock[which(stock=="NaN")] <- D
	stock
}

# One-period fleet returns
one_p_return <- function(X,S,t,p,F) {
	p[t]*S - F[t]*X
}

# Fleet pre-value function
fleet_preval <- function(X,S,D,asats,t,value_fn,p,F,igrid,...) {
	S_next <- S_(X,S,D)
	D_next <- D_(X,S,D,asats[t])
	next_state <- c(S_next,D_next)
	gridmax <- max(igrid)
	interpolation <- interpolate(next_state,igrid,value_fn)
	ifelse(next_state[2]>gridmax||L(next_state[1],next_state[2])==1,interpolation<-0,interpolation<-interpolation)
	prof <- one_p_return(X,S,t,p,F) + discount_fac*interpolation
	#if(is.infinite(prof)) {prof <- 0}
	#if(is.na(prof)) {prof <- 0}
	return(prof)
}

# Fleet pre-value function with spline interpolation
fleet_preval_spline <- function(X,S,D,asats,t,value_fn,p,F,igrid,tps_model,...) {
	S_next <- S_(X,S,D)
	D_next <- D_(X,S,D,asats[t])
	next_state <- c(S_next,D_next)
	gridmax <- max(igrid)
	interpolation <- predict(tps_model,x=cbind(S_next,D_next))
	ifelse(next_state[2]>gridmax||L(next_state[1],next_state[2])==1,interpolation<-0,interpolation<-interpolation)
	prof <- one_p_return(X,S,t,p,F) + discount_fac*interpolation
	#if(is.infinite(prof)) {prof <- 0}
	#if(is.na(prof)) {prof <- 0}
	return(prof)
}


# open access equilibrium condition
eqmcond <- function(X,S,D,fe_eqm,asats,...) {
	L(S_(X,S,D),D_(X,S,D,asats)) - fe_eqm
}

# terminal period steady state value assuming returns and costs stay constant and only replacement launches
fleet_ssval_T <- function(S,D,T,p,F,...) {
	value <- ((discount_fac^T)/(1-discount_fac))*one_p_return(L(S,D)*S,S,T,p,F)
	return(value)
} 


# # Infinite horizon satellite value
# V_ss <- function(S,D,...) {
# 	p/(1-discount_fac*(1-L(S,D)))
# }

# # Infinite horizon fleet value
# W_ss <- function(S,D,...) {
# 	(p(S)*S - F*L(S,D)*S)/(1-discount_fac)
# }
