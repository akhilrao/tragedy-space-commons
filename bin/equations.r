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

# function for "marginal survival rate" (\mathcal{L})
elL <- function(S,D,...) {
	(1 - L_prob(S,D) - S*L_S(S,D))*avg_sat_decay
}

# function for "marginal satellite return" (\alpha_1)
alpha_1 <- function(S,D,t,p,F,...) {
	ifelse(identical(F[t],numeric(0)), F_curr <- 0, F_curr <- F[t])
	ifelse(is.na(F[t]), F_curr <- 0, F_curr <- F[t])
	ifelse(is.na(p[t]), p_curr <- p[t-1], p_curr <- p[t])
	p_curr + elL(S,D)*F_curr
}

# function for "cost of marginal debris' collision" (\alpha_2)
alpha_2 <- function(S,D,t,F,...) {
	ifelse(identical(F[t],numeric(0)), F_curr <- 0, F_curr <- F[t])
	ifelse(is.na(F[t]), F_curr <- 0, F_curr <- F[t])
	0-S*L_D(S,D)*avg_sat_decay*F_curr
}

# function for "growth-launch fragment balance" (\Gamma_1)
Gamma_1 <- function(S,D,...) {
	G_S(S,D) - m*elL(S,D)
}

# function for "new fragments from current stock" (\Gamma_2)
Gamma_2 <- function(S,D,...) {
	1 - d + G_D(S,D) + m*S*L_D(S,D)
}

# fleet planner's marginal value of debris (should be exact?)
W_D <- function(S,D,t,p,F,...) {
	ifelse(identical(F[t-1],numeric(0)), F_prev <- 0, F_prev <- F[t-1])
	ifelse(identical(F[t],numeric(0)), F_prev <- 0, F_prev <- F[t-1])
	ifelse(is.na(F[t-1]), F_prev <- 0, F_prev <- F[t-1])
	ifelse(is.na(F[t]), F_prev <- 0, F_prev <- F[t-1])
	# if(t==36){
	# 		print(F_prev/discount_fac)
	# 		print(alpha_1(S,D,t,p,F))
	# 		print(Gamma_1(S,D)/Gamma_2(S,D))
	# 		print(alpha_2(S,D,t,F))
	# 		print(Gamma_1(S,D)/Gamma_2(S,D) + m)
	# 	}
	(F_prev/discount_fac - alpha_1(S,D,t,p,F) + (Gamma_1(S,D)/Gamma_2(S,D))*alpha_2(S,D,t,F))/(Gamma_1(S,D)/Gamma_2(S,D) + m)
}

# marginal external cost approximation
xi <- function(S,D,t,F,...) {
	S*L_S(S,D)*avg_sat_decay*F[t] + (Gamma_1(S,D)/Gamma_2(S,D))*alpha_2(S,D,t,F)
}

# open access equilibrium condition
eqmcond <- function(X,S,D,fe_eqm,asats,...) {
	L(S_(X,S,D),D_(X,S,D,asats)) - fe_eqm*S_(X,S,D)
	#L(S_(X,S,D),D_(X,S,D,asats)) - fe_eqm
}

# fleet planner's approximate optimality condition
optcond_approx <- function(X,S,D,fe_eqm,t,F,asats,...) {
	L_prob(S_(X,S,D),D_(X,S,D,asats[t])) - fe_eqm[t+1] + xi(S_(X,S,D),D_(X,S,D,asats[t]),t,F)/F[t]
}

# fleet planner's exact(?) optimality condition: needs d < 1 to be well-defined
optcond_exact <- function(X,S,D,clock_time,p,F,asats,...) {
	#print(W_D(S,D,clock_time,p,F))
	# print(alpha_2(S,D,clock_time,F))
	# print(Gamma_2(S,D))
	#print(S_(X,S,D))
	#print(D_(X,S,D,asats[clock_time]))
	#print(W_D(S_(X,S,D),D_(X,S,D,asats[clock_time]),(clock_time+1),p,F))
	#print(X)
	result <- W_D(S,D,clock_time,p,F) - alpha_2(S,D,clock_time,F) - discount_fac*Gamma_2(S,D)*W_D(S_(X,S,D),D_(X,S,D,asats[clock_time]),(clock_time+1),p,F)
	#print(result)
	return(result)
}

# One-period fleet returns
one_p_return <- function(X,S,t,p,F) {
	p[t]*S - F[t]*X
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

# terminal period steady state value assuming returns and costs stay constant and only replacement launches
fleet_ssval_T <- function(S,D,T,p,F,...) {
	value <- ((discount_fac^T)/(1-discount_fac))*one_p_return(L(S,D)*S,S,T,p,F)
	return(value)
} 
