#####################################################################
# Script to take a series of orbital stock levels and calibrate a physical model of debris evolution
#####################################################################
rm(list=ls())

library(ggplot2)
library(glmnet)
library(gridExtra)
library(reshape2)
library(gmm)

#dfrm <- read.csv("/home/akhil/Documents/git-repos/tragedy-space-commons/data/ST_stock_series.csv")
dfrm <- read.csv("/platter/git-repos/tragedy-space-commons/data/ST_stock_series.csv")
dfrm <- dfrm[-c(nrow(dfrm)-1,nrow(dfrm)),]

##### define functions
fitplot <- function(xvars,coefs,year,truth,title) {
	fitline <- xvars%*%coefs
	fit <- data.frame(year=year,fit=fitline,truth=truth,error=(fitline-truth))

	plot_base <- ggplot(data=fit, aes(x=year))
	plot_fitplot <- plot_base + geom_point(aes(y=truth),size=1.1) +
							geom_line(aes(y=fit),size=0.9,linetype="dashed", color="blue") +
							theme_minimal() + ggtitle(paste(title))
	plot_error <- plot_base + geom_line(aes(y=error),size=0.9) +
						geom_hline(yintercept=0,linetype="dashed") +
						theme_minimal()

	grid.arrange(plot_fitplot,plot_error,nrow=2)

}

##### Create additional X variables for regressions
dfrm$D2 <- dfrm$debris^2
dfrm$SD <- dfrm$payloads_in_orbit*dfrm$debris
dfrm$S2 <- dfrm$payloads_in_orbit^2
Sfrac <- dfrm$payloads_in_orbit/(dfrm$payloads_in_orbit + dfrm$debris)
Dfrac <- dfrm$debris/(dfrm$payloads_in_orbit + dfrm$debris)
dfrm$SSfrags <- Sfrac*dfrm$risk*dfrm$payloads_in_orbit
dfrm$SDfrags <- Dfrac*dfrm$risk*dfrm$payloads_in_orbit
#dfrm$asat_pre <- c(dfrm$num_destr_asat[-1],0)

##### Create outcome variables for regressions
D_next <- c(dfrm$debris[-1])
S_next <- c(dfrm$payloads_in_orbit[-1])
risk <- dfrm$risk[-(dim(dfrm)[1])]

##### Bind created variables into dataframe/matrix for regressions
series <- cbind(dfrm[-(dim(dfrm)[1]),],D_next,S_next,risk)
ellfit_mat <- as.matrix(subset(series,select=c(S2,SD)))
ellfit_dfrm <- subset(series,select=c(S2,SD))
sfit_mat <- as.matrix(subset(series,select=c(payloads_in_orbit,payloads_decayed,launch_successes) ))
sfit_dfrm <- subset(series,select=c(payloads_in_orbit,payloads_decayed,launch_successes) )
dfit_mat <- as.matrix(subset(series,select=c(payloads_in_orbit,payloads_decayed,debris,launch_successes,launch_failures,num_destr_asat,SSfrags,SDfrags,D2) ) )

##### Define training samples
train <- seq(from=1,to=dim(dfit_mat)[2],by=1)
#set.seed(20)
#train <- sort(floor(runif(n=15,min=1,max=58)))

##### Run regressions to calibrate things

### risk equation calibration.
risk_elnet <- glmnet(x=ellfit_mat,y=risk, alpha=0,lambda=cv.glmnet(x=ellfit_mat,y=risk,alpha=0)$lambda.min,intercept=FALSE)
coef(risk_elnet)

risk_xvars <- as.matrix(cbind(subset(series,select=colnames(ellfit_mat))))
risk_coefs <- coef(risk_elnet)[-1]
fitplot(risk_xvars,risk_coefs,series$year,risk,"Risk calibration")

risk_ols <- lm(risk ~ -1 + S2 + SD, data=ellfit_dfrm)
summary(risk_ols)

risk_coefs <- coef(risk_ols)
png(width=400,height=400,filename="../images/risk_fit_plot.png")
fitplot(risk_xvars,risk_coefs,series$year,risk,"Risk equation calibration")
dev.off()

write.csv(risk_coefs,file="calibrated_risk_eqn_coefs.csv")

### satellite equation calibration. this is a sanity check: all coefficients should be +1 or -1 to indicate accounting relationship.
sat_elnet <- glmnet(x=sfit_mat,y=S_next, alpha=1,lambda=cv.glmnet(x=sfit_mat,y=S_next,alpha=1)$lambda.min,intercept=FALSE)
coef(sat_elnet)

sat_ols <- lm(S_next ~ -1 + payloads_in_orbit + payloads_decayed + launch_successes, data=sfit_dfrm)
summary(sat_ols)

sat_parms <- coef(sat_ols)

# plot fit
sat_xvars <- as.matrix(cbind(subset(series,select=colnames(sfit_mat))))
sat_coefs <- sat_parms
fitplot(sat_xvars,sat_coefs,series$year,S_next,"Satellite accounting model")


### debris equation calibration
# Elnet fit
### The ridge fit isn't really theoretically grounded, but the coefficients all satisfy the theoretical restrictions and are plausible orders of magnitude (with the exception of delta, which seems too high). Training on about 1/4 of the sample, evenly spaced, seems to deliver a decent fit. Going with this for now as a quick-and-dirty calibration, in lieu of something more theoretically grounded from BW2009.

dfit_mat_elnet <- as.matrix(subset(dfit_mat,select=c(payloads_decayed,debris,launch_successes,num_destr_asat,SSfrags,SDfrags,D2) ) )
m2 <- glmnet(x=dfit_mat_elnet,y=D_next, alpha=0,lambda=cv.glmnet(x=dfit_mat,y=D_next,alpha=0)$lambda.min,intercept=FALSE)
coef(m2)

# GMM calibration
dfit_mat_GMM <- as.data.frame(cbind(D_next,dfit_mat_elnet))
#dfit_mat_GMM$adjusted_debris <- (dfit_mat_GMM[,3]-dfit_mat_GMM[,2])
for(i in 1:6) {dfit_mat_GMM <- cbind(dfit_mat_GMM,dplyr::lag(dfit_mat_GMM$debris,i))}
dfit_mat_GMM <- na.omit(dfit_mat_GMM)
for(i in 9:ncol(dfit_mat_GMM)) colnames(dfit_mat_GMM)[i] <- paste0("debris_l",(i-8))
#D_GMM <- dfit_mat_GMM[,1] ~ dfit_mat_GMM[,2] + dfit_mat_GMM[,3] + dfit_mat_GMM[,4] + dfit_mat_GMM[,5] + dfit_mat_GMM[,6] + dfit_mat_GMM[,7] + dfit_mat_GMM[,8]
D_GMM <- dfit_mat_GMM[,1] ~ -1 + dfit_mat_GMM[,3] + dfit_mat_GMM[,4] + dfit_mat_GMM[,5] + dfit_mat_GMM[,6] + dfit_mat_GMM[,7] + dfit_mat_GMM[,8]
res <- gmm(D_GMM,dfit_mat_GMM[,9:ncol(dfit_mat_GMM)])
res

# plot fit
m2xvars <- as.matrix(cbind(subset(series,select=colnames(dfit_mat_elnet))))
m2coefs <- coef(m2)[-1]
names(m2coefs) <- colnames(dfit_mat_elnet)
fitplot(m2xvars,m2coefs,series$year,D_next,"Debris law of motion calibration")

png(width=400,height=400,filename="../images/debris_lom_fit_plot.png")
fitplot(m2xvars,m2coefs,series$year,D_next,"Debris law of motion calibration")
dev.off()

write.csv(m2coefs,file="calibrated_debris_lom_coefs.csv")

## debris equation sensitivity analysis
## approach: draw 15/58 observations randomly 1000 times, estimate ridge on those observations, store coefficients, then calculate projections and plot
ndraws <- 1000
set.seed(20)
train_iter <- floor(runif(n= (15*ndraws),min=1,max=58))
train_iter <- matrix(train_iter,nrow=15,byrow=TRUE)
coef_iter <- matrix(0,nrow=(length(m2coefs)+1),ncol=ndraws)
rownames(coef_iter) <- c(colnames(dfit_mat_elnet),"lambda")
fitline_iter <- matrix(0,nrow=57,ncol=ndraws)

for(i in 1:ncol(train_iter)) {
	train_set <- sort(train_iter[,i])
	dfit_mat_elnet_iter <- as.matrix(subset(dfit_mat,select=c(payloads_decayed,debris,launch_successes,num_destr_asat,SSfrags,SDfrags,D2) ) )[train_set,]
	m2_iter <- glmnet(x=dfit_mat_elnet_iter,y=D_next[train_set], alpha=1,lambda=cv.glmnet(x=dfit_mat,y=D_next,alpha=1)$lambda.min,intercept=FALSE)

	coef_iter[,i] <- c(coef(m2_iter)[-1],m2_iter$lambda)
	fitline_iter[,i] <- m2xvars%*%coef_iter[-nrow(coef_iter),i]
}

mean_coefs <- rowMeans(coef_iter)

fitline_iter_aug <- as.data.frame(cbind(series$year,fitline_iter))
colnames(fitline_iter_aug)[1] <- c("year")
fitline_long <- melt(fitline_iter_aug,id='year')
names(fitline_long) <- c('year','sample','fitline')
fitline_long <- cbind(fitline_long,D_next)

debris_sensitivity <- ggplot(data=fitline_long,aes(x=year)) + 
	geom_line(aes(y=fitline,group=sample),alpha=0.015,size=0.75,color="blue") +
	geom_line(aes(y=D_next),size=1.5) +
	theme_minimal() + theme(legend.position="none") +
	geom_hline(yintercept=0,linetype="dashed") +
	ggtitle("Sensitivity of debris equation fit to sampled points") +
	xlab("year") + ylab("Debris stock") +
	geom_vline(xintercept=series$year[which(dfrm$num_destr_asat>0)],linetype="dashed")

png(width=600,height=600,filename="../images/debris_lom_sensitivity_plot.png")
debris_sensitivity
dev.off()
