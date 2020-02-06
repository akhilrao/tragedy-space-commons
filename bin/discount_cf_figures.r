##### Script to generate figures for discount rate sensitivity analysis
###
# Script flow:
# 1. Read in counterfactuals with different discount rates
# 2. Generate individual figures
# 3. Generate figure panels

#############################################################################
# 1.  Read in alternate discount rate results, make additional outcome variables
#############################################################################

files <- list.files(path="../data/counterfactuals/discount_rate/", pattern="*.csv", full.names=FALSE, recursive=FALSE)
setwd("../data/counterfactuals/discount_rate/")

container <- list()
for(i in 1:length(files)) {
	input <- read.csv(files[i])
	container[[i]] <- cbind(calc_tax_path(input),discount_rate=discount_rate_vary[i])
}

OA_OPT_full <- rbindlist(container)

#############################################################################
# 2.  Generate individual figures
#############################################################################

selection <- which(OA_OPT_full$start_time.opt==14&OA_OPT_full$discount_rate<=0.15)
OA_OPT_tax_shift <- OA_OPT_full[selection,]
shifted_tax <- OA_OPT_full[selection,]$opt_tax_path[-1]
OA_OPT_tax_shift <- OA_OPT_tax_shift[-nrow(OA_OPT_tax_shift),]
OA_OPT_tax_shift$shifted_tax <- shifted_tax

data_size <- 1

(opt_tax_path <- ggplot(data=OA_OPT_tax_shift,aes(x=year)) + 
	geom_line(aes(y=shifted_tax, group=discount_rate, color=discount_rate),size=data_size) + theme_bw() +
	scale_y_continuous(name="Optimal OUF (nominal $/sat)", labels = scales::comma) +
	scale_color_viridis(discrete=FALSE) +
	#guides(color=guide_legend(title="Planner's\ndiscount\nrate")) +
	labs(color="Planner's\ndiscount\nrate") +
	xlab("Year") +
	xlim(2020,2040) +
	ggtitle("Optimal OUF path") +
	expand_limits(y=0) +
	theme(text=element_text(family="Helvetica",size=20),
		axis.text.x=element_text(family="Helvetica",size=20),
		axis.text.y=element_text(family="Helvetica",size=20),
		plot.title=element_text(family="Helvetica",size=20),
		legend.text=element_text(family="Helvetica",size=20) ) ) 


boa_fin_dfrm.discount_cf <- OA_OPT_tax_shift[which(OA_OPT_tax_shift$year==2040),c("NPVPoA","discount_rate","year")]
row.names(boa_fin_dfrm.discount_cf) <- NULL

aes(fill=as.factor(NPVPoA)),

(boa_plot_pc.discount_cf <- ggplot(data=boa_fin_dfrm.discount_cf,aes(discount_rate,NPVPoA)) +
				geom_bar(position="dodge", stat="identity") +
				geom_hline(yintercept = 4) +
				geom_hline(,linetype="dashed",alpha=0.5,yintercept = 1) +
				ylab("Ratio of optimal NPV to BAU NPV in 2040") +
				xlab("Planner's discount rate") +
				theme_bw() +
				scale_y_continuous(minor_breaks=seq(0,40,5),breaks=c(0,1,4,10,20,30)) +
				scale_color_viridis(discrete=FALSE) +
				guides(fill=FALSE) +
				theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
				theme(text=element_text(family="Helvetica",size=20),
					axis.text.x=element_text(family="Helvetica",size=15),
					axis.text.y=element_text(family="Helvetica",size=15),
					plot.title=element_text(family="Helvetica",size=15),
					legend.text=element_text(family="Helvetica",size=15) ) )

png(width=900,height=600,filename=paste0("../../../images/discount_rate_sensitivity.png"))
	plot_grid(boa_plot_pc.discount_cf,opt_tax_path,align="v",axis="1",labels=c("a","b"),nrow=1,rel_widths=c(1/2,1/2),label_size=25)
dev.off()
