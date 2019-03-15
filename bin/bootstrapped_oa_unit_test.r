# Script to check the distribution of bootstrapped parameters and their open access simulations

rm(list=ls())

library(ggplot2)
library(gridExtra)

##### define plotting functions
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

nls_fitplot <- function(xvars,coefs,year,truth,title) {
	fitline <- xvars[,1]*(1 - exp(-coefs[1]*xvars[,1] -coefs[2]*xvars[,2]))
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

# Multiple plot function: taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Read in data
m2coefs <- read.csv("../data/calibrated_debris_lom_coefs.csv")
m2coefs_names <- as.character(m2coefs[,1])
m2coefs <- data.frame(parameters=t(c(m2coefs[,2])))
colnames(m2coefs) <- m2coefs_names

nls_coefs <- read.csv("../data/calibrated_risk_eqn_coefs.csv")[,-1]

ridge_debris_bootstrap_coefs <- read.csv("../data/bootstrapped_debris_lom_coefs.csv")
nls_risk_bootstrap_coefs <- read.csv("../data/bootstrapped_risk_eqn_coefs.csv")
B <- nrow(nls_risk_bootstrap_coefs)

## Unconditional distributions
deblom_title_hist <- ggtitle("\n")
d_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=(1-debris))) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab("Debris decay") +
		geom_vline(xintercept=(1-m2coefs[1,1]), linetype="dashed", color="blue") +
		geom_vline(xintercept=(1-mean(ridge_debris_bootstrap_coefs$debris)), linetype="dashed") +
		ggtitle(paste0("Unconditional bootstrap distributions: ", B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
		theme_bw()
m_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=launch_successes)) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab("Launch debris") +
		geom_vline(xintercept=m2coefs[1,2], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$launch_successes), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
gamma_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=num_destr_asat)) + 

		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nanti-satellite missile tests")) +
		geom_vline(xintercept=m2coefs[1,3], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$num_destr_asat), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
beta_SS_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SSfrags)) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nsatellite-satellite collisions")) +
		geom_vline(xintercept=m2coefs[1,4], linetype="dashed", color="blue") +

		geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SSfrags), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
beta_SD_hist <- ggplot(data=ridge_debris_bootstrap_coefs, aes(x=SDfrags)) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +

	geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nsatellite-debris collisions")) +
		geom_vline(xintercept=m2coefs[1,5], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(ridge_debris_bootstrap_coefs$SDfrags), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()

# dev.new()
# multiplot(deblom_title_hist, d_hist, m_hist, gamma_hist, beta_SS_hist, beta_SD_hist, 
# 		layout=matrix(c(2,1,3,4,5,6), 
# 		nrow=3, ncol=2, byrow=TRUE))


S2_hist <- ggplot(data=nls_risk_bootstrap_coefs, aes(x=S2)) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
		xlab("Satellite-satellite collision parameter") +
		geom_vline(xintercept=nls_coefs[1,1], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(nls_risk_bootstrap_coefs$S2), linetype="dashed") +
		ggtitle(paste0("Unconditional bootstrap distributions: ",B, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
		theme_bw()
SD_hist <- ggplot(data=nls_risk_bootstrap_coefs, aes(x=SD)) + 
		geom_histogram(aes(y=..density..), bins=floor(B/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
		xlab("Satellite-debris collision parameter") +
		geom_vline(xintercept=nls_coefs[1,2], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(nls_risk_bootstrap_coefs$SD), linetype="dashed") +
		ggtitle("\n") +
		theme_bw() +
		theme(panel.background = element_rect(fill = "transparent"),
			plot.background = element_rect(fill = "transparent", color = NA))

## Conditional distributions
accepted_set <- which(nls_risk_bootstrap_coefs$SD>0)
cond_ridge_debris_bootstrap_coefs <- ridge_debris_bootstrap_coefs[accepted_set,]
cond_nls_risk_bootstrap_coefs <- nls_risk_bootstrap_coefs[accepted_set,]
B_cond <- nrow(cond_ridge_debris_bootstrap_coefs)

cond_deblom_title_hist <- ggtitle("\n")
cond_d_hist <- ggplot(data=cond_ridge_debris_bootstrap_coefs, aes(x=(1-debris))) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab("Debris decay") +
		geom_vline(xintercept=(1-m2coefs[1,1]), linetype="dashed", color="blue") +
		geom_vline(xintercept=(1-mean(cond_ridge_debris_bootstrap_coefs$debris)), linetype="dashed") +
		ggtitle(paste0("Conditional bootstrap distributions: \n", B_cond, " draws (blue: original estimate, black: mean of bootstrap estimates)")) +
		theme_bw()
cond_m_hist <- ggplot(data=cond_ridge_debris_bootstrap_coefs, aes(x=launch_successes)) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab("Launch debris") +
		geom_vline(xintercept=m2coefs[1,2], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(cond_ridge_debris_bootstrap_coefs$launch_successes), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
cond_gamma_hist <- ggplot(data=cond_ridge_debris_bootstrap_coefs, aes(x=num_destr_asat)) + 

		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nanti-satellite missile tests")) +
		geom_vline(xintercept=m2coefs[1,3], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(cond_ridge_debris_bootstrap_coefs$num_destr_asat), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
cond_beta_SS_hist <- ggplot(data=cond_ridge_debris_bootstrap_coefs, aes(x=SSfrags)) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nsatellite-satellite collisions")) +
		geom_vline(xintercept=m2coefs[1,4], linetype="dashed", color="blue") +

		geom_vline(xintercept=mean(cond_ridge_debris_bootstrap_coefs$SSfrags), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()
cond_beta_SD_hist <- ggplot(data=cond_ridge_debris_bootstrap_coefs, aes(x=SDfrags)) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666",colour="dark gray") +
		xlab(paste0("Fragments from \nsatellite-debris collisions")) +
		geom_vline(xintercept=m2coefs[1,5], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(cond_ridge_debris_bootstrap_coefs$SDfrags), linetype="dashed") +
		ggtitle("\n") +
		theme_bw()

# dev.new()
# multiplot(cond_deblom_title_hist, cond_d_hist, cond_m_hist, cond_gamma_hist, cond_beta_SS_hist, cond_beta_SD_hist, 
# 		layout=matrix(c(2,1,3,4,5,6), 
# 		nrow=3, ncol=2, byrow=TRUE))

cond_S2_hist <- ggplot(data=cond_nls_risk_bootstrap_coefs, aes(x=S2)) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
		xlab("Satellite-satellite collision parameter") +
		geom_vline(xintercept=nls_coefs[1,1], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(cond_nls_risk_bootstrap_coefs$S2), linetype="dashed") +
		ggtitle(paste0("Conditional bootstrap distributions: ",B_cond, " draws\n(blue: original estimate, black: mean of bootstrap estimates)")) +
		theme_bw()
cond_SD_hist <- ggplot(data=cond_nls_risk_bootstrap_coefs, aes(x=SD)) + 
		geom_histogram(aes(y=..density..), bins=floor(B_cond/3), colour="gray", fill="white") +
		geom_density(alpha=.2, fill="#FF6666", colour="dark gray") +
		xlab("Satellite-debris collision parameter") +
		geom_vline(xintercept=nls_coefs[1,2], linetype="dashed", color="blue") +
		geom_vline(xintercept=mean(cond_nls_risk_bootstrap_coefs$SD), linetype="dashed") +
		ggtitle("\n") +
		theme_bw() +
		theme(panel.background = element_rect(fill = "transparent"),
			plot.background = element_rect(fill = "transparent", color = NA))

dev.new()
multiplot(S2_hist, SD_hist, cond_S2_hist, cond_SD_hist, 
		layout=matrix(c(1,2,3,4), 
		nrow=2, ncol=2, byrow=TRUE))

png(width=600,height=600,filename="../images/collision_rate_parameters_unconditional_and_conditional.png")
multiplot(S2_hist, SD_hist, cond_S2_hist, cond_SD_hist, 
		layout=matrix(c(1,2,3,4), 
		nrow=2, ncol=2, byrow=TRUE))
dev.off()

dev.new()
multiplot(deblom_title_hist, d_hist, m_hist, gamma_hist, beta_SS_hist, beta_SD_hist,
	cond_deblom_title_hist, cond_d_hist, cond_m_hist, cond_gamma_hist, cond_beta_SS_hist, cond_beta_SD_hist, 
		layout=matrix(c(2,3,5,1,4,6, 8,9,11,7,10,12), 
		nrow=3, ncol=4, byrow=FALSE))

png(width=1200,height=600,filename="../images/debris_lom_parameters_unconditional_and_conditional.png")
multiplot(deblom_title_hist, d_hist, m_hist, gamma_hist, beta_SS_hist, beta_SD_hist,
	cond_deblom_title_hist, cond_d_hist, cond_m_hist, cond_gamma_hist, cond_beta_SS_hist, cond_beta_SD_hist, 
		layout=matrix(c(2,3,5,1,4,6, 8,9,11,7,10,12), 
		nrow=3, ncol=4, byrow=FALSE))
dev.off()
