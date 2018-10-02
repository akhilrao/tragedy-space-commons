rm(list=ls())

library(rootSolve)
library(gridExtra)
library(ggplot2)
library(viridis)
source("equations.r")
source("simulation_functions.r")
source("simulation_algorithms.r")

#############################################################################
# Function definitions
#############################################################################



#############################################################################
# Calibration
#############################################################################

risk_cal <- read.csv("calibrated_risk_eqn_coefs.csv")
deblom_cal <- read.csv("calibrated_debris_lom_coefs.csv")
econ_series <- read.csv("econ_series.csv")

upper <- 1e15
upper_seq <- seq(from=1,by=1,length=T)
T <- 100

aSS <- risk_cal[1,2]
aSD <- risk_cal[2,2]
aDDbDD <- deblom_cal[7,2]
bSS <- deblom_cal[5,2]
bSD <- deblom_cal[6,2]
d <- deblom_cal[2,2]
Z_coef <- deblom_cal[1,2]
m <- deblom_cal[3,2]
asat_coef <- deblom_cal[4,2]

p <- 1
F <- 20 # F=60 => takes 5 years for a satellite to earn back purchase+launch cost
r_s <- p/F
r <- 0.01
discount <- 1/(1+r)
discount_fac <- discount
excess_return <- r_s-r

#############################################################################
# Generate static mb/mc pictures
#############################################################################

F <- 20 # F=60 => takes 5 years for a satellite to earn back purchase+launch cost
r_s <- p/F
r <- 0.01
#### Code to make plots for env brown bag 04/20/2018 slides

#############################################################################
# setwd and load basic functions - run this block before anything else in here
#############################################################################
discount <- 1/(1+r)
excess_return <- r_s-r

S <- 250
D <- 250

xgrid <- seq(from=0,to=1000,by=0.5)
mb <- rep(r_s,length.out=length(xgrid))
mc <- c(r+L(S_(xgrid,S,D),D_(xgrid,S,D)))
plot_mb_mc(mb,mc,xgrid)

stock_price <- rep(0.01,length.out=length(xgrid))
flow_price_t <- rep((1+r)*0.01,length.out=length(xgrid))
flow_price_t1 <- rep((1-L(S_(xgrid,S,D),D_(xgrid,S,D)))*0.01,length.out=length(xgrid))

dfrm <- data.frame(mb=mb,mc=mc,X=xgrid,stock_price=stock_price,flow_price_t=flow_price_t,flow_price_t1=flow_price_t1)

# basic plot
basic_plot <- ggplot(data=dfrm) + 
				geom_line(aes(x=X,y=mb),size=1.1) +
				geom_line(aes(x=X,y=mc),size=1.1) +
				theme_minimal() + xlab("Launches") + ylab("$")

# plot with a stock control
basic_stock_plot <- ggplot(data=dfrm, aes(x=X)) + 
						geom_line(aes(y=mb),size=1.1) +
						geom_line(aes(y=mc),size=1.1) +
						geom_line(aes(y=mc+stock_price),linetype="dashed",size=1.1) +
						theme_minimal() + xlab("Launches") + ylab("$") +
						theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

# plot with a stock control
png(file="../../images/basic_stock_plot.png",width=600,height=300)
basic_stock_plot
dev.off()

# plot with a flow control
basic_flow_plot <- ggplot(data=dfrm, aes(x=X)) + 
						geom_line(aes(y=mb),size=1.1) +
						geom_line(aes(y=mc),size=1.1) +
						geom_line(aes(y=mc+flow_price_t),linetype="dashed",size=1.1) +
						geom_line(aes(y=mb+flow_price_t1),linetype="dashed",size=1.1) +
						theme_minimal() + xlab("Launches") + ylab("$") +
						theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

# plot with a flow control
png(file="../../images/basic_flow_plot.png",width=600,height=300)
basic_flow_plot
dev.off()

#############################################################################
# Generate sample paths for intiating control
#############################################################################

aSS <- 0.00000000001
aSD <- 0.00000005
aDD <- 0.0000000001
bSS <- 100
bSD <- 75
bDD <- 1.5
d <- 0.1
m <- 1.5

F <- 15 # F=60 => takes 5 years for a satellite to earn back purchase+launch cost
r_s <- p/F
r <- 0.01
discount_fac <- 1/(1+r)
excess_return <- r_s-r
control_value <- 0.01

##### Initiate control

zeros <- rep(0,length=T)
stock_price_path <- c(rep(0,length=T/2),rep(control_value,length=T/2))
flow_price_path <- c(rep(0,length=T/2),rep(control_value,length=T/2))

stockstart <- which(stock_price_path>0)[1]-1
flowstart <- which(flow_price_path>0)[1]-1

oa_stock_ts <- oa_tsgen_stock(150,6000,T,excess_return,stock_price_path)
oa_flow_ts <- oa_tsgen_flow(150,6000,T,excess_return,flow_price_path)

stockflow <- merge(oa_stock_ts,oa_flow_ts,by=c("time"),suffix=c(".stock",".flow"))

stockflow_base <- ggplot(data=stockflow[40:70,], aes(x=time)) +
				geom_vline(xintercept=flowstart,linetype="dashed",size=1.2) +
				theme_minimal() + xlab("Time") +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- stockflow_base +
	geom_line(aes(y=launches.stock),size=1.25,color="blue") + 
	geom_line(aes(y=launches.flow),size=1.05) + 
	ylab("Launches")
satellites <- stockflow_base +
	geom_line(aes(y=satellites.stock),size=1.25,color="blue") + 
	geom_line(aes(y=satellites.flow),size=1.05) + 
	ylab("Satellites")
collision <- stockflow_base +
	geom_line(aes(y=collision_rate.stock),size=1.25,color="blue") + 
	geom_line(aes(y=collision_rate.flow),size=1.05) + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")
debris <- stockflow_base +
	geom_line(aes(y=debris.stock),size=1.25,color="blue") + 
	geom_line(aes(y=debris.flow),size=1.05) + 
	ylab("Debris")

png(file="../../images/intro_stockflow.png",width=800,height=400)
grid.arrange(launches, satellites, collision, debris, ncol=2)
dev.off()

##### Shut down orbital access

shutdown_value <- 0.06

zeros <- rep(0,length=T)
stock_price_path <- c(rep(control_value,length= (T/2 -1) ),rep(shutdown_value,length=(T/2 + 1))
flow_price_path <- c(rep(control_value,length= (T/2 -1) ),rep(shutdown_value,length=(T/2+1) ))

stockstop <- which(stock_price_path==shutdown_value)[1]-1
flowstop <- stockstop

oa_stock_ts <- oa_tsgen_stock(150,6000,T,excess_return,stock_price_path)
oa_flow_ts <- oa_tsgen_flow_shutdown(150,6000,T,excess_return,flow_price_path,shutdown_date=51)

stockflow <- merge(oa_stock_ts,oa_flow_ts,by=c("time"),suffix=c(".stock",".flow"))

stockflow_base <- ggplot(data=stockflow[40:70,], aes(x=time)) +
				geom_vline(xintercept=flowstart,linetype="dashed",size=1.2) +
				theme_minimal() + xlab("Time") +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- stockflow_base +
	geom_line(aes(y=launches.stock),size=1.25,color="blue") + 
	geom_line(aes(y=launches.flow),size=1.05) + 
	ylab("Launches")
satellites <- stockflow_base +
	geom_line(aes(y=satellites.stock),size=1.25,color="blue") + 
	geom_line(aes(y=satellites.flow),size=1.05) + 
	ylab("Satellites")
collision <- stockflow_base +
	geom_line(aes(y=collision_rate.stock),size=1.25,color="blue") + 
	geom_line(aes(y=collision_rate.flow),size=1.05) + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")
debris <- stockflow_base +
	geom_line(aes(y=debris.stock),size=1.25,color="blue") + 
	geom_line(aes(y=debris.flow),size=1.05) + 
	ylab("Debris")

png(file="../../images/shutdown_stockflow.png",width=800,height=400)
grid.arrange(launches, satellites, collision, debris, ncol=2)
dev.off()

#############################################################################
# Generate sample paths for optimal taxes
#############################################################################

fp_baseline <- fp_tsgen(10,10,100)


source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")

opt_stock_path <- optimal_stock_tax_path(10,10,100,excess_return)
oa_opt_tax <- oa_tsgen_stock(10,10,100,excess_return,opt_stock_path)



opt_flow_path <- optimal_flow_tax_path(fp_baseline$launches,10,10,excess_return)

test <- oa_tsgen_flow(10,10,100,excess_return,opt_flow_path)
test <- oa_tsgen_stock(10,10,100,excess_return,opt_stock_path)

#############################################################################
# Generate optimal stock and flow tax policy function pictures
#############################################################################

rm(list=ls())

library(plot3D)
library(viridis)

source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")
source("projection_functions.r")

rotate <- function(x) t(apply(x, 2, rev))

### economics
p <- 1
F <- 10
removal_cost <- 100
r_c <- removal_cost/F
r_s <- p/F
discount_fac <- 0.95
r <- (1 - discount_fac)/discount_fac 
fe_eqm <- r_s - r
T <- 100

### physics/engineering
d <- 0.5
m <- 3

#### statmech rate parameters - still testing
aSS <- 7.64e-6
aSD <- 1.36e-5
aDD <- 2.55e-5
bSS <- 100
bSD <- 75
bDD <- 3
sd_subs <- 2
ss_subs <- 2

nor_opt_vfn <- t(as.matrix(read.csv("nor_vfi_vfn.csv")[,-1]))
nor_opt_lpolicy <- t(as.matrix(read.csv("nor_vfi_launch_pfn.csv")[,-1]))
nor_opt_rpolicy <- t(as.matrix(read.csv("nor_vfi_removal_pfn.csv")[,-1]))

nor_oa_vfn <- t(as.matrix(read.csv("nor.oa.fleet.value.csv")[,-1]))
nor_oa_lpolicy <- t(as.matrix(read.csv("nor.oa.launch.policy.csv")[,-1]))
nor_oa_rpolicy <- t(as.matrix(read.csv("nor.oa.removal.policy.csv")[,-1]))

# set parameters based on solve
gridlist <- build_grid(0, 24, nrow(nor_oa_vfn), 1)
base_piece <- gridlist[[1]]
igrid <- gridlist[[2]]
colnames(igrid) <- c("S","D")

sat_mat <- t(matrix(igrid$S,nrow=length(base_piece),ncol=length(base_piece),byrow=FALSE))
deb_mat <- t(matrix(igrid$D,nrow=length(base_piece),ncol=length(base_piece),byrow=FALSE))

par(mfrow=c(1,2))
image2D(z=sat_mat,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal value function",cex.lab=1.5)
image2D(z=deb_mat,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal value function",cex.lab=1.5)

par(mfrow=c(1,2))
image2D(z=t(nor_opt_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal value function",cex.lab=1.5)
image2D(z=t(nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal value function",cex.lab=1.5)

opt_s_ <- (S_(nor_opt_lpolicy,sat_mat,deb_mat))
opt_d_ <- (D_(nor_opt_lpolicy,sat_mat,deb_mat))

oa_s_ <- (S_(nor_oa_lpolicy,sat_mat,deb_mat))
oa_d_ <- (D_(nor_oa_lpolicy,sat_mat,deb_mat))


par(mfrow=c(1,2))
image2D(z=t(nor_opt_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal launch policy",cex.lab=1.5)
image2D(z=t(nor_oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Open access launch policy",cex.lab=1.5)


png(file="../../images/tax_model_output.png",width=700,height=900)
par(mfrow=c(3,2))
image2D(z=t(nor_opt_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal launch rate",cex.lab=1.5)
image2D(z=t(nor_oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Open access launch rate",cex.lab=1.5)
image2D(z=opt_s_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Optimal S'[t+1])),cex.lab=1.5,contour=TRUE)
image2D(z=oa_s_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Open access S'[t+1])),cex.lab=1.5,contour=TRUE)
image2D(z=opt_d_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Optimal D'[t+1])),cex.lab=1.5,contour=TRUE)
image2D(z=oa_d_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Open access D'[t+1])),cex.lab=1.5,contour=TRUE)
dev.off()

par(mfrow=c(1,3))
image2D(z=opt_s_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Optimal S'[t+1])),cex.lab=1.5,contour=TRUE)
image2D(z=oa_s_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Open access S'[t+1])),cex.lab=1.5,contour=TRUE)
image2D(z=(opt_s_ - oa_s_),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main=expression(bold('Open access S'[t+1])),cex.lab=1.5,contour=TRUE)

write.csv(opt_s_,file="./opt_s_.csv")
write.csv(oa_s_,file="./oa_s_.csv")
write.csv((opt_s_ - oa_s_),file="opt_s_-oa_s_.csv")

opt_L_ <- L(opt_s_,opt_d_)
oa_L_ <- L(oa_s_,oa_d_)
#par(mfrow=c(1,2))

xi_ <- oa_L_*F - opt_L_*F
image2D(z=xi_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Marginal external cost of launch",cex.lab=1.5)

stock_tax <- xi_
flow_tax <- ((1+r)*0 - xi_)/(1-opt_L_)
flow_tax[is.infinite(flow_tax)] <- NA
#flow_tax[which.max(flow_tax)] <- NA
#View(flow_tax)

png(file="../../images/optimal_tax_policies.png",width=800,height=700)
par(mfrow=c(2,2))
image2D(z=opt_L_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Open access collision risk in t+1",cex.lab=1.5,cex.axis=1.5,cex.sub=1.5)
image2D(z=oa_L_,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal collision risk in t+1",cex.lab=1.5,cex.axis=1.5,cex.sub=1.5)
image2D(z=stock_tax,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal satellite tax policy",cex.lab=1.5,cex.axis=1.5,cex.sub=1.5)
image2D(z=flow_tax,x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal launch tax policy",cex.lab=1.5,cex.axis=1.5,cex.sub=1.5,colkey(font=1.5))
dev.off()

# make images of read-in policies
#png(file="../../images/det_sd_comparison.png",width=800,height=1000)
par(mfrow=c(3,2))

image2D(z=t(nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Open access fleet value",cex.lab=1.5,clab="Value",contour=TRUE)


image2D(z=t(nor_opt_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Optimal fleet value",cex.lab=1.5,contour=TRUE)

image2D(z=t(nor_opt_lpolicy - nor_oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Differences in launch policies (optimal - open access)",cex.lab=1.5)
image2D(z=t(nor_opt_vfn - nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Differences in fleet values (optimal - open access)",cex.lab=1.5)
#dev.off()



#############################################################################
# Removal calibration
#############################################################################

upper <- 1e15
upper_seq <- seq(from=1,by=1,length=T)

aSS <- 0.00000000001
aSD <- 0.00000005
aDD <- 0.0000000001
bSS <- 100
bSD <- 75
bDD <- 1.00001
d <- 0.25
m <- 1

T <- 100
p <- 1
F <- 20 # F=60 => takes 5 years for a satellite to earn back purchase+launch cost
removal_c <- 0.0015
r_s <- p/F
r_c <- removal_c/F
r <- 0.01
discount <- 1/(1+r)
excess_return <- r_s-r

#############################################################################
# Open access launching with exogenous removal
#############################################################################

rcbar_path <- c(rep(0,length=T/2),rep(0.0015,length=T/2))
Rbar_path <- c(rep(0,length=T/2),rep(10,length=T/2))
remstart <- which(Rbar_path>0)[1]-1

oa_exorem_ts_costly <- oa_tsgen_exorem(150,6000,T,excess_return,rcbar_path,Rbar_path)
oa_exorem_ts_free <- oa_tsgen_exorem(150,6000,T,excess_return,rep(0,length=T),Rbar_path)

oa_exorem_ts <- merge(oa_exorem_ts_costly,oa_exorem_ts_free,by=c("time"),suffix=c(".cost",".free"))

exorem_base <- ggplot(data=oa_exorem_ts[40:70,], aes(x=time)) +
				geom_vline(xintercept=remstart,linetype="dashed",size=1.2) +
				theme_minimal() +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- exorem_base +
	geom_line(aes(y=launches.free),size=1.25) + 
	geom_line(aes(y=launches.cost),size=1.05,color="blue") + 
	ylab("Launches") + xlab("")
satellites <- exorem_base +
	geom_line(aes(y=satellites.free),size=1.25) + 
	geom_line(aes(y=satellites.cost),size=1.05,color="blue") + 
	ylab("Satellites") + xlab("")
collision <- exorem_base +
	geom_line(aes(y=collision_rate.free),size=1.25) + 
	geom_line(aes(y=collision_rate.cost),size=1.05,color="blue") + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")  + xlab("Time")
debris <- exorem_base +
	geom_line(aes(y=debris.free),size=1.25) + 
	geom_line(aes(y=debris.cost),size=1.05,color="blue") + 
	ylab("Debris")  + xlab("Time")

png(file="../../images/intro_exorem.png",width=800,height=400)
grid.arrange(launches, satellites, collision, debris, ncol=2)
dev.off()


#############################################################################
# Open access launching with endogenous removal
#############################################################################

aSS <- 0.00000001
aSD <- 0.000005
aDD <- 0.000001
bSS <- 100
bSD <- 75
bDD <- 5

p <- 1
F <- 60 # F=60 => takes 5 years for a satellite to earn back purchase+launch cost
removal_c <- 0.0015
r_s <- p/F
r_c <- removal_c/F
r <- 0.01
discount <- 1/(1+r)
excess_return <- r_s-r

m <- 1
T <- 100
remstart <- (T/2)

D_k <- uniroot.all(kessthres,c(0.1,1000000000))

removal_c <- r_c*F
#rc_path <- c(rep(0,length=(remstart-1) ),rep(r_c,length=(T-remstart+1) ))
rc_path <- rep(r_c,length=T)
source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")
oa_endorem_ts <- oa_tsgen_endorem(1,3500,T,excess_return,rc_path,remstart,0.014)

endo_removals <- oa_endorem_ts$agg_removal[which(oa_endorem_ts$agg_removal>0)]
endo_sats <- oa_endorem_ts$satellites[which(oa_endorem_ts$agg_removal>0)]

# merge this with endorem for comparison
rcbar_path <- c(rep(0,length=T/2),rep(r_c,length=T/2))
Rbar_path <- c(rep(0,length=((T/2)-1)),endo_removals*endo_sats,0)

oa_exorem_ts_costly <- oa_tsgen_exorem(1,3500,T,excess_return,rcbar_path,Rbar_path)

rem_comp <- merge(oa_endorem_ts,oa_exorem_ts_costly,by=c("time"),suffixes=c(".endo",".exo"))

remcomp_base <- ggplot(data=rem_comp[40:70,], aes(x=time)) +
				geom_vline(xintercept=remstart,linetype="dashed",size=1.2) +
				theme_minimal() +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- remcomp_base +
	geom_line(aes(y=log(launches.exo)),size=1.25) + 
	geom_line(aes(y=log(launches.endo)),size=1.05,color="blue") + 
	ylab("ln(Launches)") + xlab("")
satellites <- remcomp_base +
	geom_line(aes(y=satellites.exo),size=1.25) + 
	geom_line(aes(y=satellites.endo),size=1.05,color="blue") + 
	ylab("Satellites") + xlab("")
collision <- remcomp_base +
	geom_line(aes(y=collision_rate.exo),size=1.25) + 
	geom_line(aes(y=collision_rate.endo),size=1.05,color="blue") + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")  + xlab("Time")
debris <- remcomp_base +
	geom_line(aes(y=debris.exo),size=1.25) + 
	geom_line(aes(y=debris.endo),size=1.05,color="blue") + 
	ylab("Debris")  + xlab("Time")

png(file="../../images/intro_endorem.png",width=800,height=400)
grid.arrange(launches, satellites, collision, debris, ncol=2)
dev.off()

endorem_base <- ggplot(data=oa_endorem_ts, aes(x=time)) +
				geom_vline(xintercept=remstart,linetype="dashed",size=1.2) +
				theme_minimal() +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- endorem_base +
	geom_line(aes(y=launches),size=1.25) + 
	ylab("Launches") + xlab("")
satellites <- endorem_base +
	geom_line(aes(y=satellites),size=1.25) + 
	ylab("Satellites") + xlab("")
collision <- endorem_base +
	geom_line(aes(y=collision_rate),size=1.25) + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")  + xlab("Time")
debris <- endorem_base +
	geom_line(aes(y=debris),size=1.25) + 
	ylab("Debris")  + xlab("Time")
removals <- endorem_base +
	geom_line(aes(y=removals),size=1.25) + 
	ylab("Aggregate removal") + xlab("")
fleetval <- endorem_base +
	geom_line(aes(y=fleet_pv),size=1.25) + 
	ylab("Fleet value") + xlab("Time")

#png(file="../../images/intro_endorem.png",width=800,height=400)
grid.arrange(launches, satellites, collision, debris, ncol=2)
#dev.off()


#############################################################################
# Optimal launching with endogenous removal
#############################################################################

# aSS <- 0.000001
# aSD <- 0.000005
# T <- 100
# r_c <- 0.00311
# removal_c <- r_c*F
#remstart <- 25
T <- 150
# rc_path <- c(rep(0,length=(remstart-1) ),rep(r_c,length=(T-remstart+1) ))
source("equations.r")
source("simulation_functions.r")
source("simulation_algorithms.r")
fp_endorem_ts <- fp_tsgen_endorem(1,3500,T,excess_return,rc_path,remstart)
#seriesgen_ts_rem(c(fp_endorem_ts$launches,fp_endorem_ts$removals,0),25,500,T)

endorem_base <- ggplot(data=fp_endorem_ts, aes(x=time)) +
				geom_vline(xintercept=remstart,linetype="dashed",size=1.2) +
				theme_minimal() +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- endorem_base +
	geom_line(aes(y=launches),size=1.25) + 
	ylab("Launches") + xlab("")
removals <- endorem_base +
	geom_line(aes(y=removals),size=1.25) + 
	ylab("Aggregate removal") + xlab("")
fleetval <- endorem_base +
	geom_line(aes(y=fleet_pv),size=1.25) + 
	ylab("Fleet value") + xlab("Time")
satellites <- endorem_base +
	geom_line(aes(y=satellites),size=1.25) + 
	ylab("Satellites") + xlab("")
collision <- endorem_base +
	geom_line(aes(y=collision_rate),size=1.25) + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")  + xlab("Time")
debris <- endorem_base +
	geom_line(aes(y=debris),size=1.25) + 
	ylab("Debris")  + xlab("Time")

#png(file="../../images/intro_endorem.png",width=800,height=400)
grid.arrange(launches, satellites, collision, removals, debris, fleetval, ncol=3)
#dev.off()

tot_endorem <- merge(oa_endorem_ts,fp_endorem_ts,by=c("time"),suffix=c(".oa",".fp"))


tot_endorem_base <- ggplot(data=tot_endorem, aes(x=time)) +
				geom_vline(xintercept=remstart,linetype="dashed",size=1.2) +
				theme_minimal() +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

launches <- tot_endorem_base +
	geom_line(aes(y=launches.oa),size=1.25) + 
	geom_line(aes(y=launches.fp),size=1.05,color="blue") + 
	ylab("Launches") + xlab("")
removals <- tot_endorem_base +
	geom_line(aes(y=agg_removal),size=1.25) + 
	geom_line(aes(y=removals),size=1.05,color="blue") + 
	ylab("Aggregate removal") + xlab("")
fleetval <- tot_endorem_base +
	geom_line(aes(y=fleet_pv.oa),size=1.25) + 
	geom_line(aes(y=fleet_pv.fp),size=1.05,color="blue") + 
	ylab("Fleet value") + xlab("Time")
satval <- tot_endorem_base +
	geom_line(aes(y=fleet_pv.oa),size=1.25) + 
	geom_line(aes(y=fleet_pv.fp),size=1.05,color="blue") + 
	ylab("Fleet value") + xlab("Time")
satellites <- tot_endorem_base +
	geom_line(aes(y=satellites.oa),size=1.25) + 
	geom_line(aes(y=satellites.fp),size=1.05,color="blue") + 
	ylab("Satellites") + xlab("")
collision <- tot_endorem_base +
	geom_line(aes(y=collision_rate.oa),size=1.25) + 
	geom_line(aes(y=collision_rate.fp),size=1.05,color="blue") + 
	geom_hline(yintercept=excess_return,linetype="dashed",color="maroon",size=0.75) +
	ylab("Expected collision rate")  + xlab("Time")
debris <- tot_endorem_base +
	geom_line(aes(y=debris.oa),size=1.25) + 
	geom_line(aes(y=debris.fp),size=1.05,color="blue") + 
	ylab("Debris")  + xlab("Time")

#png(file="../../images/optrem_v_privrem.png",width=800,height=400)
grid.arrange(launches, satellites, fleetval, removals, debris, collision, ncol=3)
#dev.off()


#############################################################################
# Removal supply and demand
#############################################################################

source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")

removal_c <- 0.0025

s_grid <- seq(from=0,to=100,length.out=15)
d_grid <- seq(from=0,to=100,length.out=100)
sd_grid <- expand.grid(s_grid,d_grid)
colnames(sd_grid) <- c("sats","debs")

startval <- 100
Loss <- L(sd_grid$sats,sd_grid$debs)
Ld <- L_D(sd_grid$sats,startval-sd_grid$debs)
MB <- Ld*s_grid*F
Ri <- rep(0,length=length(sd_grid$sats))
for(i in 1:length(sd_grid$sats)) {
 	
 	D_temp <- sd_grid$debs[i]
 	S_temp <- sd_grid$sats[i]
 	full_rem <- ifelse( S_temp>0&&D_temp>0, D_temp/S_temp, 0.01)
 	Ri[i] <- optim(par=full_rem/2, fn=satval_rem, D=D_temp, S=S_temp , launch_cost=F, removal_cost=removal_c, method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))$par
}

R <- Ri*sd_grid$sats
Dopt <- sd_grid$debs - R
#privrem(0,10,1,F,removal_c)

statmech_dfrm <- cbind(sd_grid,Loss,Ld,MB,removal_c,Ri,R,Dopt)
colnames(statmech_dfrm) <- c("sats","debs","Loss","Ld","MB","MC","Ri","R","Dopt")

statmech_mb <- ggplot(data=statmech_dfrm[which(statmech_dfrm$debs<startval),]) + 
				geom_line(aes(x=debs,y=MB,group=sats,color=sats),size=1) + 
				theme_minimal() + #theme(legend.position="none") +
				labs(title = "PMB of debris removal", x = "Debris removed", y = "Multiple of one-period return", color = "Satellites\n") +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))  +
				scale_colour_viridis()

dev.new()
statmech_mb

statmech_loss <- ggplot(data=statmech_dfrm[which(statmech_dfrm$debs<startval),]) + geom_line(aes(x=debs,y=Loss,group=sats,color=sats)) + theme_bw() + labs(title = "statmech collision rate", x = "Debris in orbit", y = "Proportion of satellites lost", color = "Satellites\n")

statmech_mbmc <- ggplot(data=statmech_dfrm[intersect(which(statmech_dfrm$debs<startval),which(statmech_dfrm$sats==100)),]) + 
					geom_line(aes(x=debs,y=MB),color="blue",size=1.2) + 
					geom_line(aes(x=debs,y=removal_c),size=1.2) + 
					theme_minimal() + 
					labs(x = "Debris removed", y = "Multiple of one-period return", color = "Satellites\n") +
					theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12)) +
				scale_colour_viridis()

statmech_mbmc

png(file="../../images/removal_mbmc.png",width=500,height=250)
statmech_mbmc
dev.off()

statmech_rid <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + geom_line(aes(x=debs,y=Ri,group=sats,color=sats)) + theme_bw() + labs(title = "Satellite owner's debris removal demand (statmech collision rate)", x = "Initial debris", y = "Removal demanded", color = "Satellites\n")
statmech_ris <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + geom_line(aes(x=sats,y=Ri,group=debs,color=debs)) + theme_bw() + labs(title = "Satellite owner's debris removal demand (statmech collision rate)", x = "Satellites", y = "Removal demanded", color = "Initial debris\n")
statmech_rd <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + 
				geom_line(aes(x=debs,y=R,group=sats,color=sats),size=1.2) + 
				theme_minimal() + 
				labs(title = "Total debris removal demanded", x = "Initial debris", y = "Removal demanded", color = "Satellites\n") +
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12)) +
				scale_colour_viridis()
statmech_rs <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + 
				geom_line(aes(x=sats,y=R,group=debs,color=debs)) + 
				theme_minimal() + 
				labs(title = "Aggregate debris removal demanded", x = "Satellites", y = "Removal demanded", color = "Initial debris\n") + 
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))
statmech_doptd <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + 
				geom_line(aes(x=debs,y=Dopt,group=sats,color=sats)) + 
				theme_minimal() + labs(title = "Cooperative post-removal debris level", x = "Initial debris level", y = "Post-removal level", color = "Number of firms\n")  +
				scale_colour_viridis()
statmech_dopts <- ggplot(data=statmech_dfrm[which(statmech_dfrm$sats>0),]) + 
				geom_line(aes(x=sats,y=Dopt,group=debs,color=debs)) + 
				theme_minimal() + labs(title = "Cooperative post-removal debris level", x = "Number of firms in orbit", y = "Post-removal level", color = "Initial debris level\n")  +
				scale_colour_viridis()

png(file="../../images/D_opt_s_d.png",width=800,height=400)
statmech_dopts
dev.off()

png(file="../../images/removal_mb_sats.png",width=800,height=400)
grid.arrange(statmech_mb,statmech_doptd,ncol=2)
dev.off()

##### Illustrate effects of changes in the cost parameters

source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")

Fgrid <- seq(from=5,to=25,by=0.1)
nsats <- 50
ndebs <- 500
startval <- 100
Loss <- L(nsats,ndebs)
Ri <- rep(-1,length=length(Fgrid))
for(i in 1:length(Fgrid)) {
 	full_rem <- ifelse( nsats>0&&ndebs>0, ndebs/nsats, 0.01)
 	Ri[i] <- optim(par=full_rem/2, fn=satval_rem, D=ndebs, S=nsats , launch_cost=Fgrid[i], removal_cost=removal_c, method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))$par
}
R <- Ri*nsats

F_statics_dfrm <- data.frame(launch_cost=Fgrid,ind_rem=Ri)

Ri_F_statics <- ggplot(data=F_statics_dfrm) + 
				geom_line(aes(x=launch_cost,y=ind_rem),size=1) + 
				theme_minimal() + 
				labs(title = "Individual cooperative debris removal demands", x = "Launch cost", y = "Removal demanded") + 
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

rcgrid <- seq(from=0.00005,to=0.006,by=0.00001)
nsats <- 50
ndebs <- 500
startval <- 100
Loss <- L(nsats,ndebs)
Ri <- rep(-1,length=length(rcgrid))
for(i in 1:length(rcgrid)) {
 	full_rem <- ifelse( nsats>0&&ndebs>0, ndebs/nsats, 0.01)
 	Ri[i] <- optim(par=full_rem/2, fn=satval_rem, D=ndebs, S=nsats , launch_cost=F, removal_cost=rcgrid[i], method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))$par
}
R <- Ri*nsats

c_statics_dfrm <- data.frame(rem_cost=rcgrid,ind_rem=Ri)

Ri_c_statics <- ggplot(data=c_statics_dfrm) + 
				geom_line(aes(x=rem_cost*F,y=ind_rem),size=1) + 
				theme_minimal() + 
				labs(title = "", x = "Removal cost", y = "") + 
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))

png(file="../../images/removal_cost_statics.png",width=800,height=400)
grid.arrange(Ri_F_statics,Ri_c_statics,ncol=2)
dev.off()

##### Illustrate the effects of different SD complementarity values

source("equations.r")

s_grid <- seq(from=1,to=100,length.out=250)
comp_grid <- seq(from=1,to=3,by=0.01)
scomp_grid <- expand.grid(s_grid,comp_grid)
colnames(scomp_grid) <- c("sats","comp")
ss_subs <- 3

ndebs <- 500
startval <- 100
Loss <- L(s_grid,ndebs)
Ri <- rep(0,length=length(scomp_grid[,1]))

for(i in 1:nrow(scomp_grid)) {
	nsats <- scomp_grid$sats[i]
	sd_subs <- scomp_grid$comp[i]
 	full_rem <- ifelse( nsats>0&&ndebs>0, ndebs/nsats, 0.01)
 	Ri[i] <- optim(par=full_rem/2, fn=satval_rem, D=ndebs, S=nsats , launch_cost=F, removal_cost=removal_c, method="L-BFGS-B", lower=0, upper=full_rem, control=list(fnscale=-1))$par
}
R <- Ri*scomp_grid$sats
scomp_dfrm <- cbind(scomp_grid,Ri,R)
scomp_high <- scomp_dfrm[which(scomp_dfrm$comp==1),]
scomp_low <- scomp_dfrm[which(scomp_dfrm$comp==2),]

Ri_sd_comp_high <- ggplot(data=scomp_high) + 
				geom_point(aes(x=sats,y=R),size=1) + 
				theme_minimal() + 
				labs(title = "High complementarity", x = "Number of firms in orbit", y = "Removal demanded") + 
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))
Ri_sd_comp_low <- ggplot(data=scomp_low) + 
				geom_point(aes(x=sats,y=R),size=1) + 
				theme_minimal() + 
				labs(title = "Low complementarity", x = "Number of firms in orbit", y = "Removal demanded") + 
				theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12))
grid.arrange(Ri_sd_comp_high,Ri_sd_comp_low,ncol=2)

Ri_sd_comp <- ggplot(data=scomp_dfrm) + 
			geom_jitter(aes(x=sats,y=R, group=comp,color=comp),size=0.75) + 
			theme_minimal() +
			labs(title = "Different degrees of complementarity", x = "Number of firms in orbit", y = "Removal demanded") + 
			theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"),plot.title=element_text(size=15),legend.title=element_text(size=12),legend.text=element_text(size=12)) +
			scale_colour_viridis()
dev.new()
Ri_sd_comp

png(file="../../images/sd_comp_vary.png",width=800,height=400)

dev.off()

##### Illustrate cost and congestion shifts

F <- 3

source("equations.r")

s_grid <- seq(from=1,to=100,length.out=250)
uk_grid <- expand.grid(s_grid,s_grid)
colnames(uk_grid) <- c("unkinds","kinds")
ndebs <- 500
Ri <- 5


L_D_ne <- function(k,S,D,Ri,...) {
	aSD*exp
}
L_DD_sm <- function(k,S,D,Ri,...) {
	0
}
L_DD_ne <- function(k,S,D,Ri,...) {
	-aDD*aDD*exp(-aSD*S-aDD*(D-k*Ri))
}
L_SD_sm <- function(k,S,D,Ri,...) {
	aSD
}
L_SD_ne <- function(k,S,D,Ri,...) {
	-aSD*aDD*exp(-aSD*S-aDD*(D-k*Ri))
}

L_D_ne <- function(k,S,D,Ri,...) {
	aSD*exp(-aSS*S -aSD*(D-k*Ri))
}

cost_shift_sm <- function(k,S,D,Ri,...) {
	-L_DD_sm(k,S,D,Ri)*k*Ri*F + L_D(S,(D-k*Ri))*F
}
cost_shift_ne <- function(k,S,D,Ri,...) {
	-L_DD_ne(k,S,D,Ri)*k*Ri*F + L_D_ne(S,(D-k*Ri))*F
}
cong_shift_sm <- function(k,S,D,Ri,...) {
	L_SD_sm(k,S,D,Ri)*k*F
}
cong_shift_ne <- function(k,S,D,Ri,...) {
	L_SD_ne(k,S,D,Ri)*k*F
}

dfrm <- data.frame(uk_grid,cost_shift_sm=cost_shift_sm(uk_grid$kinds,uk_grid$unkinds,ndebs,Ri),cong_shift_sm=cong_shift_sm(uk_grid$kinds,uk_grid$unkinds,ndebs,Ri),cost_shift_ne=cost_shift_ne(uk_grid$kinds,uk_grid$unkinds,ndebs,Ri),cong_shift_ne=cost_shift_ne(uk_grid$kinds,uk_grid$unkinds,ndebs,Ri))

cong_sm <- matrix(dfrm$cong_shift_sm,nrow=length(s_grid),ncol=length(s_grid))
cost_sm <- matrix(dfrm$cost_shift_sm,nrow=length(s_grid),ncol=length(s_grid))
cong_ne <- matrix(dfrm$cong_shift_ne,nrow=length(s_grid),ncol=length(s_grid))
cost_ne <- matrix(dfrm$cost_shift_ne,nrow=length(s_grid),ncol=length(s_grid))

image2D(z=cong_sm,x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Congestion shift")
image2D(z=cost_sm,x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Cost shift")
image2D(z=(cost_sm+cong_sm),x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Sum of cost and congestion shifts")

image2D(z=cong_ne,x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Congestion shift")
image2D(z=cost_ne,x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Cost shift")
image2D(z=(cost_ne+cong_ne),x=s_grid,y=s_grid,xlab=c("Kinds"),ylab=c("Unkinds"),col=magma(n=100),main="Sum of cost and congestion shifts")


dev.new()
Ri_sd_comp

png(file="../../images/cost_cong_shifts.png",width=800,height=400)

dev.off()


#############################################################################
# Generate value function matrix pictures and policy functions
#############################################################################
rm(list=ls())
library(ggplot2)
library(glmnet)
library(gridExtra)
library(viridis)
library(plot3D)
library(doParallel)
library(progress)

source("equations.r")
source("simulation_algorithms.r")
source("simulation_functions.r")
source("projection_functions.r")

# parameters from value function solves - need to automate this
degapprox <- 15
ncores <- 4
### economics
p <- 1
F <- 10
removal_cost <- 100
r_c <- removal_cost/F
r_s <- p/F
discount_fac <- 0.95
r <- (1 - discount_fac)/discount_fac 
fe_eqm <- r_s - r
T <- 100

### physics/engineering
d <- 0.25
m <- 3

#### statmech rate parameters - still testing
aSS <- 7.64e-6
aSD <- 1.36e-5
aDD <- 2.55e-5
bSS <- 100
bSD <- 75
bDD <- 3
sd_subs <- 2
ss_subs <- 2

##############################
##### Make SD value & policy pictures
##############################

nor_opt_vfn <- t(as.matrix(read.csv("BIGnor_vfi_vfn.csv")[,-1]))
nor_opt_lpolicy <- t(as.matrix(read.csv("BIGnor_vfi_launch_pfn.csv")[,-1]))
nor_opt_rpolicy <- t(as.matrix(read.csv("nor_vfi_removal_pfn.csv")[,-1]))

nor_oa_vfn <- t(as.matrix(read.csv("BIGnor.oa.fleet.value.csv")[,-1]))
nor_oa_lpolicy <- t(as.matrix(read.csv("BIGnor.oa.launch.policy.csv")[,-1]))
nor_oa_rpolicy <- t(as.matrix(read.csv("nor.oa.removal.policy.csv")[,-1]))

# set parameters based on solve
gridlist <- build_grid(0, 100, nrow(nor_oa_vfn), 1)
base_piece <- gridlist[[1]]
igrid <- gridlist[[2]]
colnames(igrid) <- c("S","D")
igrid_t <- translate(igrid,igrid)
igrid_basis <- as.matrix(basis(igrid_t,degapprox,isV=0))

# make images of read-in policies
png(file="../../images/det_sd_comparison.png",width=700,height=900)
par(mfrow=c(3,2))
image2D(z=t(nor_oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Open access launch policy",cex.lab=1.5,clab="Launches")
image2D(z=t(nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Open access fleet value",cex.lab=1.5,clab="Value",contour=TRUE)

image2D(z=t(nor_opt_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Optimal launch policy",cex.lab=1.5)
image2D(z=t(nor_opt_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Optimal fleet value",cex.lab=1.5,contour=TRUE)

image2D(z=t(nor_opt_lpolicy - nor_oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Differences in launch policies (optimal - open access)",cex.lab=1.5)
image2D(z=t(nor_opt_vfn - nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100), main="Differences in fleet values (optimal - open access)",cex.lab=1.5)
dev.off()

##############################
##### Make SDR value & policy pictures
##############################

opt_vfn <- t(as.matrix(read.csv("r0.15vfi_vfn.csv")[,-1]))
opt_lpolicy <- t(as.matrix(read.csv("r0.15vfi_launch_pfn.csv")[,-1]))
opt_rpolicy <- t(as.matrix(read.csv("r0.15vfi_removal_pfn.csv")[,-1]))

oa_vfn <- t(as.matrix(read.csv("oa.fleet.value.csv")[,-1]))
oa_lpolicy <- t(as.matrix(read.csv("oa.launch.policy.csv")[,-1]))
oa_rpolicy <- t(as.matrix(read.csv("oa.removal.policy.csv")[,-1]))

# set parameters based on solve
gridlist <- build_grid(0, 100, nrow(oa_vfn), 1)
base_piece <- gridlist[[1]]
igrid <- gridlist[[2]]
colnames(igrid) <- c("S","D")

sat_mat <- matrix(igrid$S,nrow=length(base_piece),ncol=length(base_piece),byrow=FALSE)
deb_mat <- matrix(igrid$D,nrow=length(base_piece),ncol=length(base_piece),byrow=FALSE)

par(mfrow=c(2,1))
image2D(z=t(deb_mat),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Cooperative removal plan",cex.lab=1.5,clab="Post-removal\ndebris level")
image2D(z=t(oa_rpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Cooperative removal plan",cex.lab=1.5,clab="Post-removal\ndebris level")

png(file="../../images/D_opt_s_d.png",width=1200,height=1100)
par(mar=c(5.5,5.5,5.5,5.5)) #c(bottom, left, top, right), clockwise from the bottom
image2D(z=( t(deb_mat)-t(oa_rpolicy)),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Cooperative removal plan",cex.lab=3.5,cex=3.5,clab="Post-removal\ndebris level")
dev.off()

# make images of read-in policies
png(file="../../images/det_sdr_comparison.png",width=700,height=600)
par(mar=c(4,4.5,2.25,3.25)) #c(bottom, left, top, right), clockwise from the bottom
par(mfrow=c(3,3))
image2D(z=t(oa_lpolicy),x=base_piece,y=base_piece,xlab=c(""),ylab=c("Satellites"),col=magma(n=100),main="Open access \nlaunch plan",cex.lab=1.5,clab="Launches")
image2D(z=t(oa_rpolicy),x=base_piece,y=base_piece,xlab=c(""),ylab=c(""),col=magma(n=100),main="Cooperative \nremoval plan",cex.lab=1.5,clab="Total \ndebris removed")
image2D(z=t(oa_vfn),x=base_piece,y=base_piece,xlab=c(""),ylab=c(""),col=magma(n=100), main="Open access-cooperative \nfleet value",cex.lab=1.5,clab="Value",contour=TRUE)

image2D(z=t(opt_lpolicy),x=base_piece,y=base_piece,xlab=c(""),ylab=c("Satellites"),col=magma(n=100),main="Optimal launch plan",cex.lab=1.5)
image2D(z=t(opt_rpolicy),x=base_piece,y=base_piece,xlab=c(""),ylab=c(""),col=magma(n=100),main="Optimal removal plan",cex.lab=1.5)
image2D(z=t(opt_vfn),x=base_piece,y=base_piece,xlab=c(""),ylab=c(""),col=magma(n=100), main="Optimal fleet value",cex.lab=1.5,contour=TRUE)

image2D(z=t(opt_lpolicy - oa_lpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Differences in launch policies \n(optimal - open access)",cex.lab=1.5)
image2D(z=t(opt_rpolicy - oa_rpolicy),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c(""),col=magma(n=100),main="Differences in removal policies \n(optimal - cooperative)",cex.lab=1.5)
image2D(z=t(opt_vfn - oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c(""),col=magma(n=100), main="Differences in fleet values \n(optimal - open access-cooperative)",cex.lab=1.5)
dev.off()






##############################################
##### Represent solved matrices as projections
##############################################

# translate to [-1,1] for chebyshev polynomials
nor_oa_vfn_t <- translate(nor_oa_vfn,nor_oa_vfn)
nor_oa_lpolicy_t <- translate(nor_oa_lpolicy,nor_oa_lpolicy)
nor_oa_rpolicy_t <- translate(nor_oa_rpolicy,nor_oa_rpolicy)
nor_opt_vfn_t <- translate(nor_opt_vfn,nor_opt_vfn)
nor_opt_lpolicy_t <- translate(nor_opt_lpolicy,nor_opt_lpolicy)
nor_opt_rpolicy_t <- translate(nor_opt_rpolicy,nor_opt_rpolicy)

# image gen code to verify translation worked well
#image2D(z=t(nor_oa_vfn),x=base_piece,y=base_piece,xlab=c("Debris"),ylab=c("Satellites"),col=magma(n=100),main="Open access fleet value")

# put matrices into panels
oa_panel <- data.frame(fleet_val=as.vector(nor_oa_vfn_t),launch_pol=as.vector(nor_oa_lpolicy),rem_pol=as.vector(nor_oa_rpolicy_t),S=igrid_t$S,D=igrid_t$D,igrid_basis)
opt_panel <- data.frame(fleet_val=as.vector(nor_opt_vfn_t),launch_pol=as.vector(nor_opt_lpolicy_t),rem_pol=as.vector(nor_opt_rpolicy_t),S=igrid_t$S,D=igrid_t$D,igrid_basis)

# create predictor and response panels
oa_launch_pol <- as.matrix(translate(subset(oa_panel,select=c(launch_pol)),igrid))
opt_launch_pol <- as.matrix(translate(subset(opt_panel,select=c(launch_pol)),igrid))

# fit models to launch policy
elnet_mix <- 0.5
oa_launch_model <- glmnet(x=igrid_basis,y=oa_launch_pol, alpha=elnet_mix,lambda=cv.glmnet(x=oa_launch_state_mat,y=oa_launch_pol,alpha=elnet_mix)$lambda.min,intercept=FALSE)
oa_launch_coefs <- as.matrix(coef(oa_launch_model)[-1])

opt_launch_model <- glmnet(x=igrid_basis,y=opt_launch_pol, alpha=elnet_mix,lambda=cv.glmnet(x=opt_launch_state_mat,y=opt_launch_pol,alpha=elnet_mix)$lambda.min,intercept=FALSE)

# generate sequential time series
oa_seq_ts <- oa_tsgen(0,0,T,fe_eqm)
oa_test_ts <- seriesgen_pfn(as.vector(nor_oa_lpolicy),0,0,igrid,T)
oa_mod_ts <- seriesgen_model(0,0,T,oa_launch_coefs,igrid,degapprox)
