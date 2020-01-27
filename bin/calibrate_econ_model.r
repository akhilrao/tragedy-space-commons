#####################################################################
# Script to take a series of satellite costs and returns data and calibrate an economic model of satellite excess returns
#####################################################################

source("plotting_functions.r")

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
aggrc$r_s_raw <- aggrc$tot_rev/aggrc$tot_cost  #calculate the raw rate of return on a satellite - the number of satellites drops out of the formula. (the number of satellites is important for this calculation since we want p and F, not pS and FS, but the Morgan Stanley report doesn't state what fleet sizes they're assuming.)
aggrc$Ft_Ft <- c(0,aggrc$tot_cost[-1]/aggrc$tot_cost[-length(aggrc$tot_cost)])
colnames(aggrc)[1] <- "year"
aggrc <- aggrc[-1,]

### scale aggregate revenues/costs by share in LEO
satellites <- read.csv("../data/ucsfaa_merged.csv")
LEO_share <- ddply(satellites, ~Launch.Year, summarize, share_in_LEO=mean(LEO), share_in_GSO=mean(GSO) )
colnames(LEO_share)[1] <- c("year")

aggrc <- merge(LEO_share,aggrc,by=c("year"))
aggrc$r_s <- aggrc$r_s_raw

dfrm <- read.csv("../data/stock_series.csv")
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]
risk <- subset(dfrm,select=c(year, risk))
risk$risk <- risk$risk/dfrm$payloads_in_orbit

dfrm <- merge(aggrc,risk,by=c("year"))
dfrm <- subset(dfrm,select=c(year,tot_rev,tot_cost,r_s,r_s_raw,risk,share_in_LEO,Ft_Ft))

### calculate implied IRR
dfrm$implied_r <- dfrm$r_s - dfrm$risk

### write out econ series dfrm
write.csv(dfrm,file="../data/econ_series.csv",row.names=FALSE)

##### model risk as a function of the return
dfrm_mat <- as.matrix(subset(dfrm,select=-c(risk,year,implied_r,tot_rev,tot_cost,r_s_raw,share_in_LEO)))
risk <- dfrm$risk

riskmodel <- lm( risk ~ r_s + Ft_Ft, data=dfrm)
summary(riskmodel)

riskparms <- coef(riskmodel)
riskxvars <- matrix(c(rep(1,length=dim(dfrm)[1]), dfrm$r_s, dfrm$Ft_Ft), ncol=3, byrow=FALSE )
# fitplot(riskxvars,riskparms,dfrm$year,dfrm$risk,title="Collision rate as a function of returns and costs","collision rate")
# dev.off()

write.csv(riskparms,file="../data/econ_series_coefs.csv")

##### calculate implied launch costs

F <- F_calc(a_1=riskparms[1],a_2=riskparms[2],a_3=riskparms[3],F_1=(dfrm$Ft_Ft[1]*dfrm$tot_cost[1]),pi_t=dfrm$tot_rev,risk=dfrm$risk)
p <- dfrm$tot_rev

printed_table_of_costs_and_returns <- data.frame(year=seq(from=aggrc$year[1],length.out=nrow(dfrm)),pi_t=dfrm$tot_rev,F_t=dfrm$tot_cost,F_hat=F)
#View(printed_table_of_costs_and_returns)

write.csv(round(printed_table_of_costs_and_returns,digits=2),file="implied_costs.csv")

##### plot objects of interest

png(width=500,height=400,filename="../images/risk_return_plot.png")
fitplot(riskxvars,riskparms,dfrm$year,dfrm$risk,title="Collision probability as a function of returns and costs","collision probability")
dev.off()
