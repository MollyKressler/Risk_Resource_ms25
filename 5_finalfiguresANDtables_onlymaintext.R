### Final Figures code for Chapter 3: Risk or Resource

## compiling code for figures in final manuscript, not including supporting information tables and figures

#####
pacman::p_load(sf, tidyverse, lubridate, ggplot2, patchwork, cowplot, flextable,MuMIn,ggplot2,lme4,stats,ggeffects,gtsummary, MCMCvis, viridis, tidybayes)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')
#####


###################
## Common dataframes and formatting
###################
    hexdata <- read.csv('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.')
    hexsf <- st_as_sf(st_read('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.shp'),crs='WGS84')%>%
      dplyr::select(-stndrd_, -stndrd_r,-sqzrisk)%>%
      rename(standard.hexshark = stndrd_hxs,
        standard.hexfish = stndrd_hxf,
        standard.hexdist2shore = stndrd_h2,
        standard.hexdistcmg = stndrd_hxd,
        standard.hexlowdensg = stndrd_hxl,
        standard.hexmeddensg = stndrd_hxm,
        standard.dist2jetty = stndrd_d2,
        standard.depth = stndrd_d,
        standard.sgPCA1 = st_PCA1,
        logit.sqzrisk = lgt_sqz,
        zlogit.sqzrisk = zlgt_sq,
        relPropPD = rlPrpPD
        )

###################
## Metadata on shark detections
###################

	meta <- read.csv('habitat_model/no ghosts of 6hr threshold/biometric_and_metadata_for_habsel_glmers_dec22.csv')%>%dplyr::select(-X)

	meta_flex<-flextable(meta%>%dplyr::select(-pit, -period)%>%
      mutate_if(is.numeric, round, digits = 3))%>%
		set_header_labels(sex = 'Sex',average_duration_hours = 'Average Duration\nbetween Detection (hours)',min_duration_hours = 'Minimum Duration\nbetween Detection (hours)',max_duration_days = 'Maximum Duration\nbetween Detection (days)',days = 'Days at\nLiberty',pcl = 'Pre-Caudal\nLength (cm)')%>%
		theme_zebra()%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
	meta_flex

	save_as_docx(meta_flex, path = 'resource_chp3/metadata_tags_inchapter3.docx')

 ## metadata table on sizes of individuals 
	m <- read_csv('lemonspredators_20192020blacktipsANDbulls/all_predators.csv', col_select=c(Acoustic.ID:Longitude))%>%
		mutate(PredatorID = seq(1,49,1))
	d <- read_csv('lemonspredators_20192020blacktipsANDbulls/predators_detectedwithinstudyarea_MKthesis2019to2020_aug2023.csv')%>%
		mutate(Acoustic.ID = as.character(elasmo))
		
	a<- left_join(d, m%>%dplyr::select(PredatorID, Acoustic.ID), relationship = 'many-to-one')%>%
		filter(!is.na(PredatorID))%>%
		group_by(PredatorID)%>%
		tally()

	m2 <- left_join(m,a)%>%	
		filter(!is.na(n))%>%
		mutate(PredatorID = seq(1,21,1)) # resetting this for table readability only
	m2	# metadata for only sharks detected within period and area of juveniles detection set 

	mean(m2$PCL)
	sd(m2$PCL) # summary stats for text blurb

	meta.flextable <- m2 %>%
		mutate('Date (ymd)' = str_split(Date,' ', simplify=TRUE)[,1])%>%
		dplyr::select(PredatorID, Species, Sex, 'Date (ymd)', PCL, Project.Location,n)%>%
		flextable()%>%
		set_header_labels(PredatorID = 'Predator ID', Project.Location = 'Original Project \n\ Location', n = 'Number \n\ of Detections')%>%
		theme_zebra()%>%
		align(align = 'center', part = 'all')%>%
		font(fontname = 'Arial', part = 'all')%>%
		fontsize(size = 10, part = 'all')%>%
		autofit()
	meta.flextable
	save_as_image(meta.flextable, 'lemonspredators_20192020blacktipsANDbulls/onlypredatorsinstudyareaandtime_predators_PCL_tagDate_tagLocation.png', webshot='webshot')

###################
## Habitat data map - from winter 2020, from Emily Courmier 
###################

	## import polygonised data from Winter 2020
	labels <- as_labeller(c(lowdensg='Low Density Seagrass',mediumdensg='Medium Density Seagrass',highdensg='High Density Seagrass',vegetated='Vegetated (terrestrial)',bare.urban='Bare &/or Urban (terrestrial)'))
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
	cmg<-st_as_sf(st_read('habitat_model/centroid_mangroves_north.shp'),crs='WSG84')
	r <- st_as_sf(st_read('receivers_in_thesis_data_dec2023.shp'),crs='WGS84')
	jettys<- st_as_sf(st_read('boatlaunchesBimini.kml'),crs='WGS84')
	nrow(jettys)

	h14 <- st_as_sf(st_read('winter2014habitat_polygonised_forSF_ECdata.shp'), crs='WGS84')
	h18 <- st_as_sf(st_read('summer2018habitat_polygonised_forSF_ECdata.shp'), crs='WGS84')
	h20 <- st_as_sf(st_read('winter2020habitat_polygonised_forSF_ECdata.shp'), crs = 'WGS84')


	## coverage of each type of habitat feature by area 
		total.area <- st_area(h20)/sum(st_area(h20))
		# values left to right: lds, mds, hds, vegetated, bare/urban
		# 0.271, 0.184, 0.127, 0.315, 0.100


	# hab maps 
	plot1 <- ggplot()+
		geom_sf(data=h14,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('2014')
	plot2 <- ggplot()+
		geom_sf(data=h18,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		theme_bw()+	
		theme(axis.text.y = element_blank(),legend.position = 'bottom', legend.direction = 'horizontal', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('2018')
	plot3 <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		theme_bw()+
		theme(axis.text.y = element_blank(),legend.position = 'bottom', legend.direction = 'horizontal', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('2020')

	plot <- (plot1 | plot2 | plot3) + plot_layout(guides = "collect") & theme(legend.position = 'bottom', text = element_text(size = 11))
	ggsave(plot,file='habitat_bimini_wint14_summ18_wint20_ECdata.png',device='png',unit='in',height=5.5,width=8.5,dpi=950)

	# hab map, same hab data as above for 2020, but with additional information on jettys, receivers etc.
	plot3b <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		geom_sf(data=r%>%filter(sharks=='no'),fill='grey82', col='white',pch=22,size=2)+
		geom_sf(data=r%>%filter(sharks=='yes'),fill='white', col='grey50',pch=21,size=2)+
		geom_sf(data=cmg,pch=23,size=3,col='#FFFFFF',fill='#032D59')+
		geom_sf(data = jettys, pch = 24, size = 1, col='goldenrod1',fill='goldenrod1')+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den. Seagrass', 'Low Den. Seagrass', 'Medium Den. Seagrass', 'Vegetated'))+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('2020')

		plot <- (plot1 | plot2 | plot3b) + plot_layout(guides = "collect") & theme(legend.position = 'bottom', legend.justification = 'centre',text = element_text(size = 11))
		# save
		ggsave(plot,file='habitat_bimini_wint14_summ18_wint20_withCoVdata4BSEM_ECdata.png',device='png',unit='in',height=5, width = 8.5, dpi=800)

	# just 2020 alone with covariates data points 
	plot3c <- ggplot()+
		geom_sf(data=h20,aes(fill=feature,col=feature),lwd=0)+
		geom_sf(data=r%>%filter(sharks=='no'),fill='grey82', col='white',pch=22,size=2)+
		geom_sf(data=r%>%filter(sharks=='yes'),fill='white', col='grey50',pch=21,size=2)+
		geom_sf(data=cmg,pch=23,size=3,col='#FFFFFF',fill='#032D59')+
		geom_sf(data = jettys, pch = 24, size = 1, col='goldenrod1',fill='goldenrod1')+
		scale_fill_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den.\n\ Seagrass', 'Low Den.\n\ Seagrass', 'Medium Den.\n\ Seagrass', 'Vegetated'))+
		scale_colour_manual(name = 'Feature',values=c('grey72','cadetblue4','cadetblue2','cadetblue3','darkolivegreen4'), labels = c('Bare/Urban', 'High Den.\n\ Seagrass', 'Low Den.\n\ Seagrass', 'Medium Den.\n\ Seagrass', 'Vegetated'))+
		theme_void()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal')+
		guides(fill=guide_legend(nrow=2,byrow=TRUE))

		ggsave(plot3c,file='solohabbimini2020_withCoVforBSEM_ECdata_themevoid.png',device='png',unit='in',width=5, height = 7,dpi=850)

###################
## BRUVS locations map
###################
		
		land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	
		data <- st_as_sf(st_read('bruvs2014AND2018_data_joinedWITHhabitat_summer18habANDwinter2014hab_july2024.shp'), crs= 'WGS84')
		### SD and HG data 
		bruvslocations <- ggplot()+
			geom_sf(data = land, col = 'grey75', lwd = 0.5)+
				geom_sf(data = data, aes(col = source, pch = source))+
				scale_color_manual(values = c('#020f75', '#7c1447'))+
				theme_bw()

		withsize <- ggplot()+
			geom_sf(data = land, col = 'grey75', lwd = 0.5)+
				geom_sf(data = data, aes(col = source, pch = source, size = maxN))+
				scale_color_manual(values = c('#020f75', '#7c1447'))+
				theme_bw()
		ggsave(withsize,file='bruvsdata_eval/bruvs_locations_hg_sd_bysource_withsize.png',device='png',unit='in',height=5, width = 5,dpi=850)


###################
## MaxN predictions
###################
	
	preds <- st_as_sf(st_read('predictions_from_simplifiedBRTs_Gerreidae_and_sppRichness_june24.shp'), crs = 'WGS84')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	maxN.plot<-ggplot()+		
		geom_sf(data=land,fill='gray75')+
		geom_sf(data=preds,aes(fill=maxN_preds),col=NA)+
		scale_fill_viridis(name='MaxN')+
		geom_sf(data=land,fill='gray75')+
		theme_bw()+
		theme(legend.position = 'bottom', legend.direction = 'horizontal')+
		theme(axis.text.y = element_blank(),legend.position = 'bottom', legend.direction = 'horizontal', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('MaxN')
	ggsave(maxN.plot,file='resource_chp3/BRTS_outputs/BRT_maxN_july2024/simplified_BRT_LOG_maxN_preds_july24.png',device='png',units='in',dpi=950,height=7,width=6)		


###################
## Descriptive plots for predator metric: relative risk at receivers
###################
  N = 35
	relp <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp'))%>%
        mutate(sqzrisk = ((relPropPD*(N-1))+0.5)/N)%>%
        mutate(logit.sqzrisk = logit(sqzrisk))%>%
        mutate(zlogit.sqzrisk = (logit.sqzrisk-mean(logit.sqzrisk))/sd(logit.sqzrisk))
    relp
   	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	
	relativePropPDdettsreceivers<-ggplot()+
	geom_sf(data=land,col='grey75',lwd=0.5)+
	geom_sf(data=relp,aes(size=ndetts,col=relPropPD))+
		scale_color_viridis(name='Relative\n\ Risk')+
		scale_size(name = 'No. of Detections')+
		theme_bw()+
		theme(legend.position = 'bottom', axis.text.x = element_text(angle = 45, vjust = 0.65))+
		guides(size = 'none')+
		ggtitle('Predation Risk')

	relativePropPDdettsreceivers
	
	ggsave(relativePropPDdettsreceivers,file='lemonspredators_20192020blacktipsANDbulls/descriptive_stats_and_figures/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.png',device='png',units='in',dpi=950,height=7,width=6)



###################
## Descriptive plots for Seagrasses PCA (raw data is in habitat map)
###################
	
	## seagrass PCA is standardised and mean-centred in hexdata
	withpca<-st_as_sf(st_read('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.shp'),crs='WGS84')
	land<-st_as_sf(st_read('bim_onlyland_noDots.kml'),crs='WGS84')	

	sgPCAplot <- ggplot()+
		geom_sf(data=land,col='grey75',lwd=0.5)+
		geom_sf(data=withpca,aes(fill=st_PCA1),col=NA)+
		scale_fill_viridis(name='Seagrass PCA')+
		theme_bw()+
		theme(legend.position = 'bottom',legend.direction = "horizontal", axis.text.x = element_text(angle = 45, vjust = 0.65))+
		ggtitle('Seagrasses')
	
	ggsave(sgPCAplot,file='resource_chp3/seagrassesPCA_distributionplot_forchpater3writeup_.png',device=png,units='in',height=8,width=4.2,dpi=850)



###################
## Three mid level predictors: seagrasses PCA, fish, relative risk 
###################

	plot<- (sgPCAplot | maxN.plot | relativePropPDdettsreceivers) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom', text = element_text(size = 8))
	ggsave(plot,file='resource_chp3/sgPCA_maxN_and_relRisk_chp3.png',device='png',units='in',dpi=950,height=4.5)		




##############################################
## Table: BRTs full vs simple ##
##############################################
	
	n1 <- readRDS('maxN_gbm_poisson_july2024.RDS') # full model
	simple.model <- readRDS('simplified_maxN_gbm_poisson_HDS_dist2shore_july2024.RDS')

	# Get n.trees and mean deviance of original & simplified models 
		mod <- list(n1, simple.model)
		
		info<-as.data.frame(matrix(ncol=3,nrow=2))
		colnames(info)=c('Model','Deviance','n.trees')
		info$Model<-c('n1', 'simple.model')
		for(i in info$Model){
			mod<-get(i)
			dev<-as.numeric(mod$cv.statistics$deviance.mean[1])
			info$Deviance[info$Model==i]<-dev
			trees<-as.numeric(mod$n.trees)
			info$n.trees[info$Model==i]<-trees
				}

		info2<-info%>%
			mutate('Explanatory variables'=c('low, medium and high density seagrass & dist. to shore','high density seagrass & dist. to shore'))%>%
				mutate_if(is.numeric, round, digits = 3)%>%
				flextable()%>%
				theme_zebra()%>%
				autofit()

		save_as_image(info2,'resource_chp3/BRTS_outputs/BRT_maxN_july2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_july2024.png',webshot='webshot')
		save_as_docx(info2, path='resource_chp3/BRTS_outputs/BRT_maxN_july2024/ntrees_and_deviance_of_full_and_simpleBRTS_chp3_july2024.docx')


##############################################
## Table: summaries of best glms ##
##############################################


# read in model RDS
	sgm <- readRDS('resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_july24.RDS')
	fm <- readRDS('resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_MaxN_July24.RDS')
	pm <- readRDS('resource_chp3/hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetrics_may24.RDS')

# make table
	models <- c(sgm, fm, pm)
	
	sgf <- flextable::as_flextable(sgm)
	pf <- flextable::as_flextable(pm)
	ff <- flextable::as_flextable(fm)
	gf <- flextable::as_flextable(gm)

	s1 <- tbl_regression(sgm,exp=FALSE,conf_level=0.95,label=list(standard.depth='Depth',standard.distcmg='Dist. to Central Mangroves',st.dist2jetty='Dist. to Jetty','st.dist2jetty:st.dist2shore'='Dist. to Jetty * Dist. to Shore'))%>%
		bold_p(t=0.05)
	f1 <- tbl_regression(fm, exp=FALSE,conf_level=0.95,label=list(st.dist2shore='Dist. to Shore',st.dist2jetty='Dist. to Nearest Jetty',standard.sgPCA1='Seagrass PCA','standard.dist2jetty:standard.hexdist2shore'='Dist. to Nearest Jetty * Dist. to Shore'))%>%
		bold_p(t=0.05)	
	p1 <- tbl_regression(pm, exp=FALSE,conf_level=0.95, label=list(standard.depth='Depth',standard.dist2jetty='Dist. to Jetty',standard.dist2shore='Dist. to Shore',standard.distcmg='Dist. to Central Mangroves', 'standard.depth:standard.dist2jetty' = ' Depth * Dist. to Jetty'))%>%
		bold_p(t=0.05)

	# side by side
	tbl_merge(tbls = list(s1,f1,p1), tab_spanner = c('Hypothesis 1', 'Hypothesis 2', 'Hypothesis 3'))
	# stacked 
	stacked <- tbl_stack(list(s1,f1,p1),group_header=c('1','2','3'))%>%
	bold_levels()
	
	show_header_names(stacked)
 
	responses <- as_tibble(x=c('Seagrasses','Fish','Large Sharks'))%>%rename(Respones=value)

	stacked.summary<-stacked%>%
		modify_header(groupname_col = '**Hypothesis**',label='**Predictor**')%>%
		as_flex_table()%>%	
		autofit()%>%
		theme_zebra()

	## AIC of three. mdoels relative to each other
	AIC(sgm, fm, pm)

	save_as_image(stacked.summary,'hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_July2024.png' ,webshot='webshot')
	save_as_docx(stacked.summary,path = 'hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_July2024.docx')
	

##############################################
## Table: pathway description, estimates, CI and Support ##
##############################################
	
	samplesList5a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter5000_burn1000_chains3_July2024.RDS')

	pathwayresults_table <- MCMCsummary(samplesList5a,round=5,pg0=TRUE,params='path', probs=c(0.05,0.95))%>%
			tibble::rownames_to_column()%>%
			rename_with(str_to_title)%>%
		rename('pg0'='P>0')%>%
	  	mutate(pg00 = case_when(Mean > 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
		arrange(-pg00)%>%
		rename(Pathway = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		mutate(Estimate = round(Estimate, 4))%>%
		mutate(lower = case_when(Path == 3 ~ round(lower,3), Path !=3 ~ round(lower,4)),upper = case_when(Path == 3 ~ round(upper,3), Path !=3 ~ round(upper,4)))%>%
		mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
		mutate('Standardised Estimate \n\ (95% credible interval)' = paste0(Estimate,' ',CI),.after='Pathway')%>%
		mutate('Path' = parse_number(Pathway), .before = 'Pathway')%>%  
		dplyr::select(-N.eff,-Rhat,-lower,-upper,-Estimate,-CI,-Sd, -pg0)%>%	
		flextable()%>%
			compose(i=1,j=2, as_paragraph('Juvenile sharks ~ Dist. to Refuge + Dist. to Shore + \n\ Seagrasses + Teleost fish'))%>%
			compose(i=2,j=2, as_paragraph('Juvenile sharks ~ Depth + Dist. to Shore +  \n\ Dist. to Jetty + Dist. to Refuge\n\ + Depth*Dist. to Jetty + Relative Risk'))%>%
			compose(i=3,j=2, as_paragraph('Juvenile sharks ~ Dist. to Jetty + Dist. to Shore + \n\ Dist. to Refuge + Seagrasses'))%>%
			align(j=3:4, align = 'center', part = 'all')%>%
			align(j=1, align = 'center', part = 'all')%>%
			font(fontname = 'Arial', part = 'all')%>%
			color(color='black',part='all')%>%
			fontsize(size = 10, part = 'all')%>%
			autofit()%>%
			theme_zebra()
	pathwayresults_table

	save_as_image(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model5a_niter5000_burn1000_chains3_July2024.png')	
	save_as_docx(pathwayresults_table,path='resource_chp3/nimblemodel_outputs/pathwayresultssummary_model5a_niter5000_burn1000_chains3_July2024.docx')	

##############################################
## Caterpillars: path estimates BSEM ##
##############################################
    samplesList6a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.RDS')

     # grab draws with gather_draws and create label for paths based on iterations and sequence of paths minN to maxN. 
    
    d6 <- gather_draws(samplesList6a,path[])%>%
      group_by(.chain)%>%
      mutate(pathID = paste0('path',rep(1:3, each=8000)))%>% 
      mutate(pathIDnum = rep(1:3, each=8000))%>%
      ungroup() # each path estimate for each chain (of 3) in order, starting with path[1] first estimate in chain 1 
    
    # use tidybayes to plot 
    
    caterpillars <- ggplot(d6, aes(y=reorder(pathID,pathIDnum,decreasing=T),x=.value))+
      stat_pointinterval(.width=c(.50,.8),point_size=2)+
      ylab('Path ID')+
      xlab('Parameter Estimate & CI')+
      geom_vline(xintercept=0,linetype=3)+
      guides(col = 'none')+
      theme_bw()
    caterpillars  
    ggsave(caterpillars,file='resource_chp3/nimblemodel_outputs/caterpillarsPlot__model6_niter5000_burn1000_chains3_jan2025.png',device='png',dpi=400,width=5,height=4,units='in')
    

##############################################
## Plots from predictions of effects from path 3 into 2020 hexdata ##
##############################################
	
	## load data 
  hexsf <- st_as_sf(st_read('data_for_bayes_structural_EQ_modelling_DF2_HEXAGONpredictions_andBRTs_july24.shp'),crs='WGS84')%>%
    rename(standard.hexshark = stndrd_hxs,
      standard.hexfish = stndrd_hxf,
      standard.hexdist2shore = stndrd_h2,
      standard.hexdistcmg = stndrd_hxd,
      standard.hexlowdensg = stndrd_hxl,
      standard.hexmeddensg = stndrd_hxm,
      standard.dist2jetty = stndrd_d2,
      standard.depth = stndrd_d,
      standard.sgPCA1 = st_PCA1,
      zlogit.sqzrisk = zlgt_sq,
      relPropPD = rlPrpPD
      )

	land <- st_as_sf(st_read('bim_onlyland_noDots.kml'), crs = 'WGS84')

	p3pred <- readRDS('resource_chp3/path_inference/path3_means_andHDI_at_hexagons_model5a_july2024.RData')
	
	## join spatial geomtry by jcode to preds
	
		p3.sf <- st_as_sf(left_join(hexsf, p3pred, by='jcode'))%>%
		  dplyr::select(jcode, Mean, Low, Upp, geometry)

	## path 3 plots - Absolute values for mean, upper, lower effect estimates
	p3.lower <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = (Low)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
	  theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(-1, 0.5), oob = scales::squish, name = 'Prediction')+
	  theme(axis.text.x = element_text(angle=45, hjust = 1),  legend.direction = "horizontal", legend.position = "bottom")+
    labs(subtitle = 'Lower')
	 p3.lower

	p3.upper <- ggplot()+
		geom_sf(data = p3.sf, aes(fill = (Upp)), lwd=0)+
		geom_sf(data = land, col = 'grey75', lwd=0.5)+
	  theme_bw()+
    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(-1, 0.5), oob = scales::squish, name = 'Prediction')+
	  theme(axis.text.x = element_text(angle=45, hjust = 1), legend.direction = "horizontal", legend.position = "bottom")+
    labs(subtitle = 'Upper')
    p3.upper

  # upper and lower side by side 
    p3upper.noyaxis <- p3.upper + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())
    p3.hdis.wide.axescollected <- (p3.lower | p3upper.noyaxis) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
 		ggsave(p3.hdis.wide.axescollected, file = 'resource_chp3/path_inference/path3_HDI_WIDEpanel__estimates_spatial_july24.png', device = 'png', unit = 'in', dpi = 900, width = 8)

	## path 3 prediction plot, mean, with 2 receivers of highest predaiton risk score and the jetties marked 

		relp <- st_as_sf(st_read('lemonspredators_20192020blacktipsANDbulls/relativepredatorrisk_at_receivers_April2019December2020_lemonsANDblacktips.shp'))%>%
			slice_max(relPropPD, n = 2)
		jettys<- st_as_sf(st_read('boatlaunchesBimini.kml'),crs='WGS84')


	  p3.mean.withlocations <- ggplot()+
	    geom_sf(data = p3.sf, aes(fill = (Mean)), lwd=0)+
	    geom_sf(data = land, col = 'grey75', lwd=0.5)+
	    geom_sf(data = relp, col = 'purple3', pch = 19, size = 4.25)+
			geom_sf(data = jettys, pch = 24, size = 3, col='goldenrod1',fill='goldenrod1')+
	    theme_bw()+
	    scale_fill_gradientn(colors = c('#fafcfc', '#073B46'), limits = c(-1,0.5), oob = scales::squish, name = 'Mean')+
	      theme(axis.text.x = element_text(angle=45, hjust = 1))
	   ggsave(p3.mean.withlocations, file = 'resource_chp3/path_inference/path3_mean_estimates_spatial_withtop2RelRisk_jetties_model5_july24.png', device = 'png', unit = 'in', dpi = 900, width = 5)
	 



















##