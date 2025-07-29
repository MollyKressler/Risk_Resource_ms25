## Simplified Boosted Regression Tree Model for Teleost Families MaxN compared to Full, and predictions into spatially referenced dataframes

## created by Molly M Kressler 


###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd/')

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,arm,webshot2,sfdep, patchwork, cowplot, tidybayes, stats)

#####
## - Load data 
#####

	# shapefiles, useful for plotting
		land<-st_as_sf(st_union(st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')))
		
		# data post August 2023 - use this
		hab.grid<-st_as_sf(st_read('winter2020habitat_hexagon_grid_ECdata.shp'),crs='WGS84')%>%mutate(jcode=as_factor(jcode))
		hab.noland<-st_as_sf(st_read('winter2020habitat_hexagon_grid_noLAND_ECdata.shp'),crs='WGS84')%>%
			mutate(jcode=as_factor(jcode))

			
	# Data frames and shape files with data 
		## need a df and a sf for predicting into: jcodes and habitat data, and Season - W and D
		sf4preds<-st_as_sf(st_read('winter2020habitat_hexagon_grid_NOland_ECdata.shp'),crs='WGS84')
		b <-read.csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')


	# Model RDS - BRTS, with all 9 predictor variables, 'full'
		n1 <- readRDS('resource_chp3/model_RDS/maxN_gbm_poisson_oct2024.RDS')

		# deprecated june 2024 & july 2024
			sppfull<-readRDS('resource_chp3/model_RDS/SW_Species_brt_tc3_lr001_gaussian_aug23.RDS')
			famfull<-readRDS('resource_chp3/model_RDS/SW_Families_brt_tc9_lr001_gaussian_aug23.RDS')
			gerrfull<-readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_aug23.RDS')		
			richfull <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc5_lr0001_poisson_NoSeason_jun24.RDS')
			richfull_tc3 <- readRDS('resource_chp3/model_RDS/Species_Rchness_Poisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')
			gerrfull <- readRDS('resource_chp3/model_RDS/gerrPoisson_brt_tc3_lr0001_poisson_NoSeason_jun24.RDS')



######
## - Simplify the Full Model 
######
	# Simplify - 

	fullmodel<-n1

	simple<-gbm.simplify(fullmodel,n.drops=3) # 4 predictors, must have at least 2, so n.drops = 2

	simple.model<-gbm.step(b,gbm.x=simple$pred.list[[1]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 
	simple1.model<-gbm.step(b,gbm.x=simple$pred.list[[2]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 
	simple2.model<-gbm.step(b,gbm.x=simple$pred.list[[3]],gbm.y=fullmodel$gbm.call$gbm.y,tree.complexity=fullmodel$gbm.call$tree.complexity,learning.rate=as.numeric(fullmodel$gbm.call$learning.rate),bag.fraction=fullmodel$gbm.call$bag.fraction,family=fullmodel$gbm.call$family,plot.main = TRUE) 

	simple$pred.list
	simple.model$cv.statistics  # dist2shore and hds
	simple1.model$cv.statistics  # dist2shore and hds # best model without tide and low density seagrass. We want to keep tide
	simple2.model$cv.statistics  # dist2shore and hds

	saveRDS(simple1.model,'resource_chp3/model_RDS/simplified_maxN_gbm_poisson_MDS_HDS_dist2shore_oct2024.RDS')

	# Add tide back in because knowledge suggests it is important

	simple.withtide <- gbm.step(n,gbm.x=c('tide','prp_mds','prp_hds','dist2shore'),gbm.y='maxN',tree.complexity=3,learning.rate=0.001,bag.fraction=0.75,family='poisson',plot.main = TRUE) 
	saveRDS(simple.withtide,'resource_chp3/model_RDS/simple_maxN_gbm_poisson_TIDE_MDS_HDS_dist2shore_oct2024.RDS')

	# Get n.trees and mean deviance of original & simplified models 
		mod <- list(n1, simple1.model, simple.withtide)
		
		info<-as.data.frame(matrix(ncol=3,nrow=3))
		colnames(info)=c('Model','Deviance','n.trees')
		info$Model<-c('n1', 'simple1.model', 'simple.withtide')
		for(i in info$Model){
			mod<-get(i)
			dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
			info$Deviance[info$Model==i]<-dev
			trees<-as.numeric(mod$n.trees)
			info$n.trees[info$Model==i]<-trees
				}

		info2<-info%>%
			mutate('Explanatory variables'=c('tide state, low, medium and high density seagrass & dist. to shore','high density seagrass & dist. to shore', 'tide state, medium and high density seagrass & dist. to shore'), 'Model' = c('Full', 'Simplified', 'Simple with Tide'))%>%
				mutate_if(is.numeric, round, digits = 3)%>%
				flextable()%>%
				theme_zebra()%>%
				autofit()

		save_as_image(info2,'resource_chp3/BRTS_outputs/BRT_maxN_oct2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_oct2024.png',webshot='webshot')
		save_as_docx(info2, path='resource_chp3/BRTS_outputs/BRT_maxN_oct2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_oct2024.docx')



######
## - Evaluate simplified model 
######
		simplified <- readRDS('resource_chp3/model_RDS/simplified_maxN_gbm_poisson_MDS_HDS_dist2shore_oct2024.RDS')
		simple.model <- readRDS('resource_chp3/model_RDS/simple_maxN_gbm_poisson_TIDE_MDS_HDS_dist2shore_oct2024.RDS') # including tide based on knowledge 
		hab.labels<-(c('dist2shore'='Dist. to Shore (m)', 'prp_lds'='Prop. of \n\ Low Density \n\ Seagrass','prp_mds'='Prop. of  \n\ Medium Density \n\ Seagrass','prp_hds'='Prop. of  \n\ High Density \n\ Seagrass', 'tide' = 'Tide State'))

		infl.simple<-simple.model$contributions
		simple.relinf<-infl.simple%>%
			mutate(var=fct_reorder(var,rel.inf))%>%
			ggplot(aes(x=var,y=rel.inf,fill=rel.inf))+
				geom_bar(stat='identity')+
				scale_fill_distiller(direction=1,palette='Blues',limits=c(1,48),guide='none')+
				theme_bw()+
				coord_flip()+
				scale_x_discrete(labels=hab.labels)+
				ylab('Relative Influence')+
				xlab(NULL)	+
				theme(text = element_text(size = 14))
		ggsave(simple.relinf,file='resource_chp3/BRTS_outputs/BRT_maxN_oct2024/relative_influence_vars_maxN_BRT_SIMPLE_withTide_oct2024.png',device='png',units='in',height=4,width=5.5,dpi=900)
		
		res_simple <- as_tibble(resid(simple.model))%>%rename(resids = value)
		F1 <- as_tibble(predict(simple.model))%>%rename(fitted = value)
		diag_simple <- bind_cols(F1, res_simple)

		resVfit_simple <- ggplot(data = diag_simple,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()
		res_simple_hist <- ggplot(data = res_simple, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(res_n1)/100,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_simple_qq <- ggplot(data=F1, aes(sample = fitted))+
			stat_qq(size=1,pch=21,col = '#3A6C74', fill = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_simple <- resVfit_simple+res_simple_hist+res_simple_qq
		
		ggsave(diagnostics_simple, file = 'resource_chp3/BRTS_outputs/BRT_maxN_oct2024/diagnostics_BRT_SIMPLE_withTide_maxN_oct2024.png', device = 'png', unit = 'in', height = 4, width = 8, dpi = 850)	

######
## - Predict into 2020 habitat data - hexagon grid
######
	simple.model <- readRDS('resource_chp3/model_RDS/simple_maxN_gbm_poisson_TIDE_MDS_HDS_dist2shore_oct2024.RDS') # including tide based on knowledge 

	sf4preds_L <- sf4preds %>% mutate(tide = 'L')
	sf4preds_H <- sf4preds %>% mutate(tide = 'H')
	sf4preds <- bind_rows(sf4preds_L,sf4preds_H)%>%
		mutate(tide = as_factor(tide))
	
	maxN.preds<-predict.gbm(simple.model,sf4preds,n.trees=simple.model$gbm.call$best.trees,type='response')
	
	maxN.preds<-as.data.frame(maxN.preds) 

	predicts.sf<-bind_cols(sf4preds, maxN.preds)

	#####
	## - save
	#####
		st_write(predicts.sf,'predictions_hex_from_simplifiedBRTs_maxN_oct24.shp',driver='ESRI Shapefile',append=FALSE) # append command set to FALSE replaces the old file with the new one. If set to TRUE it would add the new data as further layers.
		st_write(predicts.sf,'predictions_hex_from_simplifiedBRTs_maxN_oct24.csv', driver = 'CSV')

	## plot

	preds <- st_as_sf(st_read('predictions_hex_from_simplifiedBRTs_maxN_oct24.shp'), crs = 'WGS84')
		
	maxN.plotLOG<-ggplot()+geom_sf(data=preds,aes(fill=log(maxN_preds)),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,5),guide=guide_colourbar(title=' MaxN of 4\n\ Families'))+
		geom_sf(data=land,fill='gray98')+
		facet_wrap(~tide)+
		theme_bw()+
		theme(strip.background = element_blank(), strip.placement = "outside")+
		scale_y_continuous(limits = c(25.66, 25.78), breaks = c(25.68,25.72,25.76))+
		scale_x_continuous(limits = c(-79.31,-79.24), breaks = c(-79.3,-79.27, -79.24))
	ggsave(maxN.plotLOG,file='resource_chp3/BRTS_outputs/simplified_BRT_LOG_maxN_preds_oct24.png',device='png',units='in',height=4.5,width=6,dpi=850)		

	maxN.plot<-ggplot()+geom_sf(data=preds,aes(fill=maxN_preds),col=NA)+
		scale_fill_distiller(palette='YlOrRd',direction=1,limits=c(0,100),guide=guide_colourbar(title=' MaxN of 4\n\ Families'))+
		geom_sf(data=land,fill='gray98')+
		facet_wrap(~tide)+
		theme_bw()+
		theme(strip.background = element_blank(), strip.placement = "outside")+
		scale_y_continuous(limits = c(25.66, 25.78), breaks = c(25.68,25.72,25.76))+
		scale_x_continuous(limits = c(-79.31,-79.24), breaks = c(-79.3,-79.27, -79.24))
	ggsave(maxN.plot,file='resource_chp3/BRTS_outputs/simplified_BRT_maxN_preds_oct24.png',device='png',units='in',height=4.5,width=6,dpi=850)

######
## - Predict into 2020 habitat data - buffers
######

	buffs.sf <- st_as_sf(st_read('receivers_in_thesis_data_dec2023_withEC2020Habdata.shp'),crs='WGS84')%>%
		rename(prp_hds = prop_hdsg)
	bl <- buffs.sf %>% mutate(tide = 'L')
	bh <- buffs.sf %>% mutate(tide = 'H')
	buffs.sf <- bind_rows(bl, bh) %>% 
		mutate(tide = as_factor(tide))%>%
		rename(prp_mds = prop_mdsg)

	simple.model <- readRDS('resource_chp3/model_RDS/simple_maxN_gbm_poisson_TIDE_MDS_HDS_dist2shore_oct2024.RDS')

	maxN.preds<-predict.gbm(simple.model,buffs.sf,n.trees=simple.model$gbm.call$best.trees,type='response')
	maxN.preds<-as.data.frame(maxN.preds)

	predicts.sf<-bind_cols(buffs.sf, maxN.preds)

	#####
	## - save
	##### 
		st_write(predicts.sf,'sf_with_predictions_INTO_BUFFERS_fromBRTs_oct2024.shp',driver='ESRI Shapefile',append=FALSE) 

		st_write(predicts.sf,'sf_with_predictions_INTO_BUFFERS_fromBRTs_oct2024.csv', driver = 'CSV')


######
## - Compare raw data to predictions 
######

	pred <- st_as_sf(st_read('predictions_hex_from_simplifiedBRTs_maxN_oct24.shp'), crs = 'WGS84')%>%
		mutate(jcode = as.character(jcode))
	b <-read_csv('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.csv')%>%
		dplyr::select(jcode, maxN)%>%
		mutate(jcode = as.character(jcode))

	# join 
		data <- left_join(pred, b, by = 'jcode', relationship = 'many-to-many')

	# plot
	boxplot<- ggplot()+
		geom_boxplot(data = pred, aes(y = log(maxN_preds+1)), fill = '#3A6C74', alpha = 0.2)+
		geom_jitter(data = b, aes(x = 0, y = log(1+maxN)), pch = 19, alpha = 0.75,,col = '#3A6C74')+
		theme_bw()+
		theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
		labs(y = 'MaxN (log+1)', x = NULL)+
		ggtitle('Simplified Predictions versus Observed')+
		facet_wrap(~tide)+
		theme(strip.background = element_blank(), strip.placement = "outside")
	ggsave(boxplot,file='resource_chp3/BRTS_outputs/BRT_maxN_oct2024/predsVSraw_maxN_simplifiedBRT_oct24.png',device='png',units='in',height=4,width=4.5,dpi=850)






