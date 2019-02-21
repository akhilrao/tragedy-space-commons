#####################################################################
# Script to take a series of satellite costs and returns data and calibrate an economic model of satellite excess returns
#####################################################################

rm(list=ls())

library(ggplot2)
library(glmnet)
library(plyr)
library(gridExtra)

# function to plot estimated objects
fitplot <- function(xvars,coefs,year,truth,title) {
	fitline <- xvars%*%coefs
	colnames(fitline)=NULL
	fit <- data.frame(year=year,fit=fitline,truth=truth,error=(fitline-truth))

	plot_base <- ggplot(data=fit, aes(x=year))
	plot_fitplot <- plot_base + geom_line(aes(y=truth),size=1.1) +
							geom_line(aes(y=fit),size=0.9,linetype="dashed", color="blue") +
							theme_minimal() + ggtitle(paste(title))
	plot_error <- plot_base + geom_line(aes(y=error),size=0.9) +
						geom_hline(yintercept=0,linetype="dashed") +
						theme_minimal()

	grid.arrange(plot_fitplot,plot_error,nrow=2)

}

# function to calculate sequence of economic launch costs recursively
F_calc <- function(a_1,a_2,a_3,F_1,pi_t,risk,...) {
	F <- rep(0,length=length(risk))
	F[1] <- F_1
	for(i in 1:(length(risk)-1)) {
		F[i+1] <- (a_2*pi_t[i] + a_3*F[i])/(risk[i+1] - a_1)
	}
	return(F)
}

##### read in data, aggregate, and reshape 
aggrc <- read.csv("../data/industry_revenues.csv")
aggrc <- subset(aggrc, select=c(Category, Commercial.Infrastructure.and.Support.Industries,Ground.Stations.and.Equipment, Satellite.Manufacturing..Commercial., Satellite.Launch.Industry..Commercial., Insurance.Premiums, Commercial.Space.Products.and.Services)) 
aggrc <- ddply(aggrc, ~Category, transform, tot_rev = Commercial.Space.Products.and.Services, tot_cost = (Commercial.Infrastructure.and.Support.Industries + Ground.Stations.and.Equipment + Satellite.Manufacturing..Commercial. + Satellite.Launch.Industry..Commercial.) )
aggrc$r_s_raw <- aggrc$tot_rev/aggrc$tot_cost  #calculate the raw rate of return on a satellite - the number of satellites drops out of the formula 
aggrc$Ft_Ft <- c(0,aggrc$tot_cost[-1]/aggrc$tot_cost[-length(aggrc$tot_cost)])
colnames(aggrc)[1] <- "year"
aggrc <- aggrc[-1,]

### scale aggregate revenues/costs by share in LEO
satellites <- read.csv("../data/ucsfaa_merged.csv")
LEO_share <- ddply(satellites, ~Launch.Year, summarize, share_in_LEO=mean(LEO), share_in_GSO=mean(GSO) )
colnames(LEO_share)[1] <- c("year")

aggrc <- merge(LEO_share,aggrc,by=c("year"))
#aggrc$r_s <- aggrc$r_s_raw*aggrc$share_in_LEO # disabled because (1) I don't have the right measure of share; (2) it feels like data snooping to use the share of launches instead of share of satellites; (3) shouldn't OLS "learn" the share and embed it in the parameter estimates?
aggrc$r_s <- aggrc$r_s_raw

dfrm <- read.csv("../data/ST_stock_series.csv")
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]
risk <- subset(dfrm,select=c(year, risk))
risk$risk <- risk$risk/dfrm$payloads_in_orbit

dfrm <- merge(aggrc,risk,by=c("year"))
dfrm <- subset(dfrm,select=c(year,tot_rev,tot_cost,r_s,r_s_raw,risk,share_in_LEO,Ft_Ft))

### calculate implied IRR
dfrm$implied_r <- dfrm$r_s - dfrm$risk

### write out econ series dfrm
write.csv(dfrm,file="econ_series.csv",row.names=FALSE)

##### model risk as a function of the return
dfrm_mat <- as.matrix(subset(dfrm,select=-c(risk,year,implied_r,tot_rev,tot_cost,r_s_raw,share_in_LEO)))
risk <- dfrm$risk

riskmodel <- lm( risk ~ r_s + Ft_Ft, data=dfrm)
summary(riskmodel)

riskparms <- coef(riskmodel)
riskxvars <- matrix(c(rep(1,length=dim(dfrm)[1]), dfrm$r_s, dfrm$Ft_Ft), ncol=3, byrow=FALSE )
fitplot(riskxvars,riskparms,dfrm$year,dfrm$risk,title="Collision rate as a function of returns and costs")
dev.off()

write.csv(riskparms,file="econ_series_coefs.csv")

##### calculate implied launch costs

F <- F_calc(a_1=riskparms[1],a_2=riskparms[2],a_3=riskparms[3],F_1=(dfrm$Ft_Ft[1]*dfrm$tot_cost[1]),pi_t=dfrm$tot_rev,risk=dfrm$risk)
p <- dfrm$tot_rev

printed_table_of_costs_and_returns <- data.frame(year=seq(from=aggrc$year[1],length.out=nrow(dfrm)),pi_t=dfrm$tot_rev,F_t=dfrm$tot_cost,F_hat=F)
#View(printed_table_of_costs_and_returns)

write.csv(round(printed_table_of_costs_and_returns,digits=2),file="implied_costs.csv")

##### plot objects of interest

png(width=400,height=400,filename="../images/risk_return_plot.png")
fitplot(riskxvars,riskparms,dfrm$year,dfrm$risk,title="Collision probability as a function of returns and costs")
dev.off()
















# # use NLS to estimate implied path of stock constraint shadow price
# launch_path <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/ST_stock_series.csv")

# sp_NLS <- as.data.frame(launch_path)
# sp_NLS <- merge(dfrm,sp_NLS,by=c("year","risk"))

# SP_NLS <- function(parms,pi_t,F_t,risk) {
# 	T <- length(risk)
# 	momcon <- rep(1,length=T)
# 	shadow_price <- parms[2:(T+1)]
# 	momcon[1] <- risk[1] - pi_t[1]/(F_t[1]) + parms[1]- shadow_price[t]
# 	for(t in 1:(T-1)){
# 		momcon[t+1] <- risk[t+1] - pi_t[t+1]/(F_t[t+1]) + parms[1]*F_t[t]/(F_t[t+1]) - shadow_price[t+1]
# 	}
# 	# print("Shadow prices are: ")
# 	# print(paste0(shadow_price))

# 	#print(paste0("Discount rate is: ", format( (parms[1]-1),digits=3)))
# 	parms[2:(T+1)] <- shadow_price
# 	objective <- sum(momcon^2) #t(momcon)%*%momcon
# 	#print(paste0("Objective value is: ", format(objective,digits=7)))
# 	return(objective)
# }

# parms <- c(1,rep(1,length=length(sp_NLS$risk)))
# #SP_NLS(parms,pi_t=sp_NLS$tot_rev,F_t=sp_NLS$tot_cost,risk=sp_NLS$risk)
# sp_res <- optim(par=parms,fn=SP_NLS,lower=0,pi_t=sp_NLS$tot_rev,F_t=sp_NLS$tot_cost,risk=sp_NLS$risk,method="L-BFGS-B")

# shadow_price <- sp_res$par[-1]
# imp_avg_r <- sp_res$par[1] - 1
# print("Shadow prices are: ")
# print(paste0(shadow_price))
# print(paste0("Discount rate is: ", format( imp_avg_r,digits=3) ))
# print(paste0("Objective value is: ", format(sp_res$value,digits=7)))
	
# sp_NLS <- cbind(sp_NLS,shadow_price)

# launches <- ggplot(sp_NLS,aes(x=year)) + geom_line(aes(y=(launch_successes+launch_failures)),size=1.1) +
# geom_line(aes(y=launch_successes),color="blue",size=1.1) +
# 							theme_minimal() + ggtitle("Total launches (attempts/year)")
# shadow_prices <- ggplot(sp_NLS,aes(x=year)) + geom_line(aes(y=shadow_price),size=1.1) +
# 							theme_minimal() + ggtitle("Price of implied stock control ($100m/year)")
# revs_costs <- ggplot(sp_NLS,aes(x=year)) + geom_line(aes(y=tot_rev),size=1.1,color="blue") +
# geom_line(aes(y=tot_cost),size=1.1,color="red") +
# 							theme_minimal() + ggtitle("Revenues and costs ($100m/year)")
# risk <- ggplot(sp_NLS,aes(x=year)) + geom_line(aes(y=risk),size=1.1) +
# 							theme_minimal() + ggtitle("Risk (inst. prob.)")
# grid.arrange(launches,revs_costs,risk,shadow_prices,ncol=1)



# # use NLS to estimate implied path of launch constraint shadow price
# launch_path <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/ST_stock_series.csv")

# lp_NLS <- as.data.frame(launch_path)
# lp_NLS <- merge(dfrm,lp_NLS,by=c("year","risk"))

# P_NLS <- function(parms,pi_t,F_t,risk) {
# 	T <- length(risk)
# 	momcon <- rep(1,length=T)
# 	shadow_price <- parms[2:(T+1)]
# 	#momcon[1] <- risk - (1 - F_t[1]/F_t[2]+shadow_price[2])) - pi_t[2]/(F_t[2]+shadow_price[2]) + parms[1]*F_t[1]/(F_t[2]+shadow_price[2])
# 	for(t in 1:(T-1)){
# 		#momcon[t+1] <- risk[t+1] - (1 - shadow_price[t]/(F_t[t+1]+shadow_price[t+1])) - pi_t[t+1]/(F_t[t+1]+shadow_price[t+1]) + parms[1]*F_t[t]/(F_t[t+1]+shadow_price[t+1])
# 		momcon[t+1] <- risk[t+1] - (1/(F_t[t+1])) - pi_t[t+1]/(F_t[t+1]) + parms[1]*F_t[t]/(F_t[t+1])
# 	}
# 	# print("Shadow prices are: ")
# 	# print(paste0(shadow_price))

# 	#print(paste0("Discount rate is: ", format( (parms[1]-1),digits=3)))
# 	parms[2:(T+1)] <- shadow_price
# 	objective <- sum(momcon^2) #t(momcon)%*%momcon
# 	#print(paste0("Objective value is: ", format(objective,digits=7)))
# 	return(objective)
# }

# parms <- c(1,rep(1,length=length(lp_NLS$risk)))
# #P_NLS(parms,pi_t=lp_NLS$tot_rev,F_t=lp_NLS$tot_cost,risk=lp_NLS$risk)
# res <- optim(par=parms,fn=P_NLS,lower=0,pi_t=lp_NLS$tot_rev,F_t=lp_NLS$tot_cost,risk=lp_NLS$risk,method="L-BFGS-B")

# shadow_price <- res$par[-1]
# imp_avg_r <- res$par[1] - 1
# print("Shadow prices are: ")
# print(paste0(shadow_price))
# print(paste0("Discount rate is: ", format( imp_avg_r,digits=3) ))
# print(paste0("Objective value is: ", format(res$value,digits=7)))
	
# lp_NLS <- cbind(lp_NLS,shadow_price)

# launches <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=(launch_successes+launch_failures)),size=1.1) +
# geom_line(aes(y=launch_successes),color="blue",size=1.1) +
# 							theme_minimal() + ggtitle("Total launches (attempts/year)")
# shadow_prices <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=shadow_price),size=1.1) +
# 							theme_minimal() + ggtitle("Implied price of launch constraint ($100m/year)")
# revs_costs <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=tot_rev),size=1.1,color="blue") +
# geom_line(aes(y=tot_cost),size=1.1,color="red") +
# 							theme_minimal() + ggtitle("Revenues and costs ($100m/year)")
# risk <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=risk),size=1.1) +
# 							theme_minimal() + ggtitle("Risk (inst. prob.)")
# grid.arrange(launches,revs_costs,risk,shadow_prices,ncol=1)


# ### random initial condition draws

# shadow_price_dist <- matrix(0,nrow=300,ncol=length(lp_NLS$risk))
# imp_avg_r_dist <- rep(0,length.out=300)
# objective_dist <- rep(0,length.out=300)

# for(i in 1:300) {
# 	#parms <- rep(runif(min=0,max=100,n=1),length=(length(lp_NLS$risk)+1) )
# 	parms <- runif(min=0,max=1,n=(length(lp_NLS$risk)+1) )
# 	res <- optim(par=parms,fn=P_NLS,lower=0,pi_t=lp_NLS$tot_rev,F_t=lp_NLS$tot_cost,risk=lp_NLS$risk,method="L-BFGS-B")
# 	shadow_price_dist[i,] <- res$par[-1]
# 	imp_avg_r_dist[i] <- res$par[1] - 1
# 	objective_dist[i] <- res$value
# }

# distns <- cbind(shadow_price_dist,imp_avg_r_dist,objective_dist)

# summary(distns)

# best_parms <- distns[which.min(distns[,12]),]
# best_sp <- best_parms[1:length(lp_NLS$risk)]

# best_fit <- cbind(lp_NLS,best_sp)

# launches <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=(launch_successes+launch_failures)),size=1.1) +
# geom_line(aes(y=launch_successes),color="blue",size=1.1) +
# 							theme_minimal() + ggtitle("Total launches (attempts/year)")
# best_shadow_prices <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=best_sp),size=1.1,color="blue") +
# 							geom_line(aes(y=shadow_price),size=1.1) +
# 							theme_minimal() + ggtitle("Implied price of launch constraint ($100m/year)")
# revs_costs <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=tot_rev),size=1.1,color="blue") +
# geom_line(aes(y=tot_cost),size=1.1,color="red") +
# 							theme_minimal() + ggtitle("Revenues and costs ($100m/year)")
# risk <- ggplot(lp_NLS,aes(x=year)) + geom_line(aes(y=risk),size=1.1) +
# 							theme_minimal() + ggtitle("Risk (inst. prob.)")
# grid.arrange(launches,revs_costs,risk,best_shadow_prices,ncol=1)



