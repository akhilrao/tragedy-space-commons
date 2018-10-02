#####################################################################
# Script to take a series of satellite costs and returns data and calibrate an economic model of satellite excess returns
#####################################################################

rm(list=ls())

library(ggplot2)
library(glmnet)
library(plyr)
library(gridExtra)


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


##### read in data, aggregate, and reshape 
aggrc <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/industry_revenues.csv")
aggrc <- subset(aggrc, select=c(Category, Commercial.Infrastructure.and.Support.Industries,Ground.Stations.and.Equipment, Satellite.Manufacturing..Commercial., Satellite.Launch.Industry..Commercial., Insurance.Premiums, Commercial.Space.Products.and.Services)) 
aggrc <- ddply(aggrc, ~Category, transform, tot_rev = Commercial.Space.Products.and.Services, tot_cost = (Commercial.Infrastructure.and.Support.Industries + Ground.Stations.and.Equipment + Satellite.Manufacturing..Commercial. + Satellite.Launch.Industry..Commercial.) )
aggrc$r_s_raw <- aggrc$tot_rev/aggrc$tot_cost
aggrc$r_c <- c(0,aggrc$tot_cost[-1]/aggrc$tot_cost[-length(aggrc$tot_cost)])
colnames(aggrc)[1] <- "year"
aggrc <- aggrc[-1,]

### scale aggregate revenues/costs by share in LEO
satellites <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/ucsfaa_merged.csv")
LEO_share <- ddply(satellites, ~Launch.Year, summarize, share_in_LEO=mean(LEO), share_in_GSO=mean(GSO) )
colnames(LEO_share)[1] <- c("year")

aggrc <- merge(LEO_share,aggrc,by=c("year"))
aggrc$r_s <- aggrc$r_s_raw*aggrc$share_in_LEO

dfrm <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/ST_stock_series.csv")
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]
risk <- subset(dfrm,select=c(year, risk))

dfrm <- merge(aggrc,risk,by=c("year"))
dfrm <- subset(dfrm,select=c(year,tot_rev,tot_cost,r_s,r_s_raw,risk,share_in_LEO,r_c))

### calculate implied IRR
dfrm$implied_r <- dfrm$r_s - dfrm$risk

### write out econ series dfrm
write.csv(dfrm,file="econ_series.csv",row.names=FALSE)

##### model risk as a function of the return
dfrm_mat <- as.matrix(subset(dfrm,select=-c(risk,year,implied_r,tot_rev,tot_cost,r_s_raw,share_in_LEO)))
risk <- dfrm$risk
riskelnet <- glmnet(x=dfrm_mat,y=risk,alpha=0,lambda=cv.glmnet(x=dfrm_mat,y=risk,alpha=0)$lambda.min,intercept=TRUE)
riskelnet
coef(riskelnet)

riskmodel <- lm( risk ~ r_s + r_c, data=dfrm)
summary(riskmodel)

riskparms <- coef(riskmodel)
riskxvars <- matrix(c(rep(1,length=dim(dfrm)[1]), dfrm$r_s, dfrm$r_c), ncol=3, byrow=FALSE )

##### plot objects of interest

png(width=400,height=400,filename="../images/risk_return_plot.png")
fitplot(riskxvars,riskparms,dfrm$year,dfrm$risk,title="ECOB risk index as a function of returns and costs")
dev.off()
