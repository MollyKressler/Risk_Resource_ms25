## Boosted Regression Tree Model for Teleost Families MaxN, full model

## created by Molly M Kressler 


###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd/')

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,arm,webshot2,sfdep,sp,spdep,beepr, HDInterval, patchwork, cowplot, tidybayes)


###################################
##########     END     ############
###################################

############ ############ ############ ############ ############ 
# Boosted Regression Trees 
############ ############ ############ ############ ############ 

	n <-read.csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')%>%
		mutate(tide = as_factor(tide))
	head(n)

######
## - RUN MODELS 
######
## use gbm.step() to both cross validate the number of trees to use, and then run that model. 

	n1 <-gbm.step(n,gbm.x=c('tide','prp_lds','prp_mds','prp_hds','dist2shore'),gbm.y='maxN',tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='poisson',plot.main = TRUE)  # 2800 tress

	saveRDS(n1,'resource_chp3/model_RDS/maxN_gbm_poisson_oct2024.RDS')

	## output of gbm.step explained 
		# $fitted = the fitted values form the final tree. on the response scale. 
		# $fitted.vars = the variance of the fitted values, on the resposne scale
		# $residuals = the residuals of the fitted values, on the response scale. 
		# $contributons  = relative importance of the variables. 
		# $self.statistics = evaluation statistics, calculated on fitted values. Shows 'fit' of the model on the training data. NOT to be reported as model performance. 
		# $cv.statistics = most appropriate for evaluation. 
		# $weights = the weights used in fitting the models. defualt is 1. 
		# cv.values = the mean of CV (cross validating) estmates of predictive deviance. Calculated at each stagewise step. This is used in the auto-generated plots of tress versus deviance. 
		# $cv.loss.ses = standard errors in CV estimates of predictive deviance at each step in the stagewise process
		# $cv.roc.matrix = matrix of the values for area under the curve estimated on the excluded data, instead of deviance in the cv.loss.matrix.


######
## - EVALUATE MODELS
######
	# read in model objects 
		n1 <- readRDS('resource_chp3/model_RDS/maxN_gbm_poisson_oct2024.RDS')

	# extract the data on relative influence of each variable. and plot it. 
		hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass'))

		infl.n<-n1$contributions
		n.relinf<-infl.n%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Blues',limits=c(0.1,40),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)+
				theme(text = element_text(size = 14))
		ggsave(n.relinf,file='resource_chp3/BRTS_outputs/BRT_maxN_oct2024/relative_influence_vars_maxN_oct2024.png',device='png',units='in',height=4,width=5.5,dpi=900)
	
	# diagnostics plots
		res_n1 <- as_tibble(resid(n1))%>%rename(resids = value)
		F1 <- as_tibble(predict(n1))%>%rename(fitted = value)
		diag_n1 <- bind_cols(F1, res_n1)

		resVfit_n1 <- ggplot(data = diag_n1,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()
		res_n1_hist <- ggplot(data = res_n1, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(res_n1)/100,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_n1_qq <- ggplot(data=F1, aes(sample = fitted))+
			stat_qq(size=1,pch=21,col = '#3A6C74', fill = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_n1 <- resVfit_n1+res_n1_hist+res_n1_qq
		
		ggsave(diagnostics_n1, file = 'resource_chp3/BRTS_outputs/BRT_maxN_oct2024/diagnostics_BRT_maxN_oct2024.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	



