##### functions to be used in deterministic satellite-debris model dynamic programming solver script

# Collision probability 
L <- function(S,D,...) {
	#1 - exp(-aSS*S-aSD*D)			#negexp rate
	#SD <- ifelse(S*D==0,1e-7,S*D)
	pmax(pmin(aSS*S^ss_subs +aSD*(S*D)^sd_subs,1),0)		#statmech rate
}

# Collision probability D derivative
L_D <- function(S,D,...) {
	ifelse(L(S,D)<1, aSD*S, 1)		#statmech rate
}

# Collision probability S derivative
L_S <- function(S,D,...) {
	ifelse(L(S,D)<1, 2*aSS*S + aSD*D, 1)		#statmech rate
}

# debris growth function
G <- function(S,D,...) {
	sat_caused <- ifelse( (S+D)>0, L(S,D)*(S/(S+D)), 0)
	deb_caused <- ifelse( (S+D)>0, L(S,D)*(D/(S+D)), 0)
	#newfrags <- bSS*sat_caused*S + bSD*deb_caused*S + bDD*(1 - exp(-aDD*D))*D 				#negexp rate
	newfrags <- bSS*sat_caused*S + bSD*deb_caused*S + bDD*aDD*D^2	#statmech rate
	return(newfrags)
}

# Satellite law of motion 
S_ <- function(X,S,D,...) {
	(1 - L(S,D))*S + X
}

# Debris law of motion
D_ <- function(X,S,D,...) {
	stock <- D*(1-d) + G(S,D) + m*X
	stock[which(stock=="NaN")] <- D
	stock
}

# Kessler threshold
kessthres <- function(D,...) {
	G(0,D) - d*D
}

# Fleet returns
one_p_return <- function(X,S,...) {
	p*S - F*X
}

# Fleet returns with removal
one_p_return_rem <- function(X,R,S,...) {
	p*S - removal_cost*R - F*X
}

# Infinite horizon satellite value
V_ss <- function(S,D,...) {
	p/(1-discount_fac*(1-L(S,D)))
}

# Infinite horizon fleet value
W_ss <- function(S,D,...) {
	(p*S - F*L(S,D)*S)/(1-discount_fac)
	#(p + (1 - L(S))*F)*S/(1-discount_fac)
}

# open access equilibrium condition
eqmcond <- function(X,S,D,fe_eqm,...) {
	L(S_(X,S,D),D_(X,S,D)) - fe_eqm
}

# open access equilibrium condition with exogenous removal
eqmcond_exorem <- function(X,S,D,fe_eqm,r_cbar,Rbar,...) {
	L(S_(X,S,D),D_(X,S,D)) - fe_eqm + r_cbar*Rbar
}

# open access equilibrium condition with a stock control
eqmcond_stock <- function(X,S,D,fe_eqm,stock,...) {
	fe_eqm - L(S_(X,S,D),D_(X,S,D)) + stock
}

# open access equilibrium condition with a flow control
eqmcond_flow <- function(X,S,D,fe_eqm,flow_t,flow_t1,...) {
	fe_eqm - L(S_(X,S,D),D_(X,S,D)) - (1+r)*flow_t + (1 - L(S_(X,S,D),D_(X,S,D)))*flow_t1
}

# open access satellite value with removal - maximize this directly (instead of using the FOC) to avoid issues with multiple optima, UNLESS those can be ruled out (should be fine for statmech?)
satval_rem <- function(Ri,S,D,launch_cost,removal_cost,...) {
	value <- p - removal_cost*Ri + (1-L(S,D-S*Ri))*launch_cost
	value <- ifelse(is.na(value),p - removal_cost*(Ri-0.00005) + (1-L(S,D-S*(Ri-0.00005) ))*launch_cost,value)
	#print(c(Ri,value))
	return(value)
}

satval_rem2 <- function(R,S,D,launch_cost,removal_cost,...) {
	value <- p - removal_cost*R/S + (1-L(S,D-R))*launch_cost
	value <- ifelse(is.na(value),p - removal_cost*((R/S)-0.00005) + (1-L(S,D-(R-0.00005) ))*launch_cost,value)
	#print(c(Ri,value))
	return(value)
}

remcond <- function(S,D,removal_cost,launch_cost,...) {
	removal_cost - L_D(S,D)*S*launch_cost
}

# marginal external cost calculation: difference in loss rates
mec <- function(oa_xsd,fp_xsd) {
	oa_loss <- L(S_(oa_xsd[,1],oa_xsd[,2],oa_xsd[,3]),D_(oa_xsd[,1],oa_xsd[,2],oa_xsd[,3]))
	fp_loss <- L(S_(fp_xsd[,1],fp_xsd[,2],fp_xsd[,3]),D_(fp_xsd[,1],fp_xsd[,2],fp_xsd[,3]))
	mec <- oa_loss - fp_loss
	return(mec)
}

# terminal period steady state value
fleet_ssval_T <- function(S,D,T,...) {
	#value <- ((discount_fac^T)/(1-discount_fac*(1-L(S,D))))*one_p_return(L(S,D)*S,S)
	value <- ((discount_fac^T)/(1-discount_fac))*one_p_return(L(S,D)*S,S)
	return(value)
} 

# terminal period steady state value with removal
fleet_ssval_T_rem <- function(R,S,D,T,...) {
	#value <- ((discount_fac^T)/(1-discount_fac*(1-L(S,D))))*one_p_return(L(S,D)*S,S)
	value <- ((discount_fac^T)/(1-discount_fac))*one_p_return_rem(L(S,D)*S,R,S)
	return(value)
} 