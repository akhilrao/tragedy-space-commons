##### functions to be used in deterministic satellite-debris model dynamic programming solver script

# Collision rate 
L <- function(S,D,...) {
	S*(1-exp( -(aSS*S) -(aSD*D) ))
	#pmax(pmin(aSS*S^2 + aSD*(S*D),S),0)		#1 statmech rate -- treat as total number of collisions rather than P(arbitrary single collision)
}

# Collision probability
L_prob <- function(S,D,...) {
	(1-exp( -(aSS*S) -(aSD*D) ))
	#pmax(pmin(aSS*S^2 + aSD*(S*D),S),0)		#1 statmech rate -- treat as total number of collisions rather than P(arbitrary single collision)
}

# Derivatives of collision rate
L_S <- function(S,D,...) {
	aSS*exp(- aSS*S - aSD*D)
	#ifelse(L_prob(S,D)==1, 0, aS + 2*aSS*S + aSD*D)
}

L_D <- function(S,D,...) {
	aSD*exp(- aSS*S - aSD*D)
	#ifelse(L_prob(S,D)==1, 0, aD + 2*aDD*D + aSD*S)
}

# debris growth function
G <- function(S,D,...) {
	#aSS*bSS*S^2 + aSD*bSD*S*D + aDDbDD*D^2
	bSS*(1-exp(-aSS*S))*S + bSD*(1-exp(-aSD*D))*S
}

# average new fragments function S derivative
G_S <- function(S,D,...) {
	bSS*(1-exp(-aSS*S)) + bSS*aSS*S*exp(-aSS*S) + bSD*(1-exp(-aSD*D))
}

# average new fragments function D derivative
G_D <- function(S,D,...) {
	aSD*bSD*S*exp(-aSD*S)
}

# function for "marginal survival rate" (\mathcal{L})
elL <- function(S,D,...) {
	(1 - L_prob(S,D) - S*L_S(S,D))*avg_sat_decay
}

# function for "marginal satellite return" (\alpha_1)
alpha_1 <- function(S,D,t,...) {
	p[t] + elL(S,D)*F[t]
}

# function for "cost of marginal debris' collision" (\alpha_2)
alpha_2 <- function(S,D,t,...) {
	-S*L_D(S,D)*avg_sat_decay*F[t]
}

# function for "growth-launch fragment balance" (\Gamma_1)
Gamma_1 <- function(S,D,...) {
	G_S(S,D) - m*elL(S,D)
}

# function for "new fragments from current stock" (\Gamma_2)
Gamma_2 <- function(S,D,...) {
	1 - d + G_D(S,D) + m*S*L_D(S,D)
}

# function for the user cost
xi <- function(S,D,t,...) {
	S*L_S(S,D)*avg_sat_decay*F + (Gamma_1(S,D)/Gamma_2(S,D))*alpha_2(S,D,t)
}

# Satellite law of motion 
S_ <- function(X,S,D,...) {
	#S*avg_sat_decay*(1-L(S,D)) + X
	(S - L(S,D))*avg_sat_decay + X #treating L(S,D) as the number of satellites lost in collisions, rather than probability single satellite is lost. assume decay happens after collisions.
}

# Debris law of motion
D_ <- function(X,S,D,asats,...) {
	#stock <- D*(1-d) + G(S,D) + m*X + Z_coef*avg_sat_decay*S + asat_coef*asats
	stock <- D*(1-d) + G(S,D) + m*X + asat_coef*asats
	stock[which(stock=="NaN")] <- D
	stock[which(stock=="NA")] <- D
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
	return(prof)
}

# open access equilibrium condition
eqmcond <- function(X,S,D,fe_eqm,asats,...) {
	L(S_(X,S,D),D_(X,S,D,asats)) - fe_eqm*S_(X,S,D)
	#L(S_(X,S,D),D_(X,S,D,asats)) - fe_eqm
}

# fleet planner's approximate optimality condition
optcond <- function(X,S,D,fe_eqm,t,...) {
	L_prob(D_(X,S,D),S_(X,S,D)) - fe_eqm + xi(D_(X,S,D),S_(X,S,D),t)/F[t]
}

# terminal period steady state value assuming returns and costs stay constant and only replacement launches
fleet_ssval_T <- function(S,D,T,p,F,...) {
	value <- ((discount_fac^T)/(1-discount_fac))*one_p_return(L(S,D)*S,S,T,p,F)
	return(value)
} 
