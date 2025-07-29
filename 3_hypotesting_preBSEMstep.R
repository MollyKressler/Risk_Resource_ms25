## Testing associations between ecological nodes

## inspired by methods in Bennett et al. Ecology & Evolution

## created by Molly M Kressler, updated Feb 2025 for model6

# Three hypotheses tested with the 'dredge' function in the package 'MuMIn' (Multi-Model inference)
https://rdrr.io/cran/MuMIn/man/dredge.html 

#############################
## - Load Workspace and Data
#############################

pacman::p_load(tidyverse,MuMIn,ggplot2,flextable,cowplot,patchwork,lme4,stats,ggeffects,gtsummary)
  options(na.action = "na.fail")

setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')

pointdata<-read.csv('pointdata_juvlemons_withAllCoV_MKthesis20192020.csv')%>%
	mutate(st_tide = case_when(tidephs == 'L' ~ 0, tidephs ==	'H' ~ 1)) 

hexdata<-read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%
		mutate(jcode=as.numeric(jcode))%>%
	mutate(st_tide = case_when(tidephs == 'L' ~ 0, tidephs ==	'H' ~ 1)) 

#############################
## - Define Hypotheses
#############################

	##Define some hypotheses to test:
	## Determine which ecological nodes make up the lowest hierarchy of the path analysis, and feed into the mid-level (seagrasses, prey distirbution, and predator pressure)

	#	1.  Seagrasses and distance metrics: st_SG ~ depth x distcmg x dist2jetty x dist2shore x tide

	#	2.  Fish, tide (L/H), distance abiotics, and medium and high density seagrasses:  fish ~ depth * dist2jetty * dist2shore * distcmg * med seagrass * high seagrass * tide

	#	3. Large Shark detections, tide, abiotics:  global model: predation ~ depth * shore * refuge * jettys * tide 


#############################
## - Hypothesis 1: seagrasses PCA
#############################

	sg <- glm(st_SG ~ st_depth * st_d2sh * st_dstc * st_d2jetty , data=hexdata) #define the global model

	sg.dredge <- dredge(sg, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model
	
	sg.dredged.models <- get.models(sg.dredge,subset=TRUE)
	summary(get.models(sg.dredge,1)[[1]]) #"best" model
	best.sg <- get.models(sg.dredge,1)[[1]]

	sg.modelsranked.tabled <- sg.dredge %>%
	  as_tibble %>%
	  round(3)%>%
	  mutate(weight=round(weight,3),
	         model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2, 3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  save_as_image(sg.modelsranked.tabled,'resource_chp3/hypotesting_dredge_results/dredged_results_seagrass_response_depth_dist2shore_distcmg_dist2jetty_Oct2024.png',webshot='webshot')


#############################
## - Hypothesis 2 : Fish and abiotics 
#############################

	hexdata <- hexdata %>% mutate(as.factor(tidephs))

	fish.glm <- glm(st_maxN ~ st_depth * st_d2sh * st_dstc * st_d2jetty + st_mds + st_hds + st_tide, data=hexdata) #define the global model

	fish.dredge <- dredge(fish.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,5)) # dredge from the global model

	fish.modelsranked.tabled <- fish.dredge%>%
		as_tibble %>%
		round(3)%>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	 filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()

	  save_as_image(fish.modelsranked.tabled,'resource_chp3/hypotesting_dredge_results/dredged_results_fishHEXresponse_Oct2024.png',webshot='webshot')
	get.models(fish.dredge,1)[[1]]


#############################
## - Hypothesis 3:  Large sharks
#############################

  ## predators, point
	pred.glm <- glm(st_risk ~ st_tide * st_depth * st_d2shore * st_distcmg + st_d2jetty, data=pointdata) #define the global model
	
	pred.dredged <- dredge(pred.glm, beta='sd', evaluate=TRUE, trace=FALSE, extra='R^2',m.lim=c(0,4)) # dredge from the global model

	press.glm.modelsranked.tabled <- pred.dredged%>%
		as_tibble %>%
		round(3) %>%
	  mutate(weight=round(weight,9), model = 1:n()) %>%
	  mutate(Null = ifelse(df == 2, 1, NA))%>%
	  pivot_longer(cols=starts_with(c('st','Null')),values_to='estimate')%>%	  
	  filter(!is.na(estimate)) %>%
	  group_by(pick(2,3,4,5,6,7,8)) %>%
	  summarise(model = paste(name, collapse = ' + ')) %>%
	  ungroup() %>%
	  dplyr::select(model, everything())%>% 
	  arrange(-weight) %>%
	  flextable()%>%	  
	  theme_zebra()%>%
	  set_header_labels(model = 'Model',delta = 'dAICc')%>%
	  align(align = 'left', part = 'all')%>%
	  color(color='black',part='all')%>%
	  fontsize(size = 10, part = 'all')%>%
	  autofit()
	  
	  save_as_image(press.glm.modelsranked.tabled,'resource_chp3/hypotesting_dredge_results/dredged_results_RelPropPDresponse_depth_dist2shore_distcmg_dist2jetty_tide_oct2024.png',webshot='webshot')

#############################
## - GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# H1: seagrasses and distance metrics

	sgm <- get.models(sg.dredge,1)[[1]]

# H2: fish and distance metrics
		fm <- get.models(fish.dredge,1)[[1]]
		fm
	# deprecated
		fm 	<- glm(standard.hexfish ~ standard.hexdist2jetty + standard.hexdistcmg + standard.sgPCA + standard.hexdist2jetty*standard.hexdistcmg, data=hexdata, family=gaussian)

		gm <- glm(standard.hexgerr ~ standard.hexdist2jetty + standard.hexdistcmg + standard.sgPCA + standard.hexdist2jetty*standard.hexdistcmg, data=hexdata, family=gaussian) # SD gerries only 

# H3: large sharks and distance metrics

	pm <- 	get.models(pred.dredged,1)[[1]]

# save the models as RDS 
	saveRDS(sgm,'resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_oct2024.RDS')
	saveRDS(fm,'resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_MaxN_oct2024.RDS')
	saveRDS(pm,'resource_chp3/hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetricsoct2024.RDS')

#############################
## - Summary table of GLM/GLMMs of 'best' models: the effect of distance metrics on predictor variables
#############################

# read in model RDS
	sgm <- readRDS('resource_chp3/hypotesting_dredge_results/seagrasses_glm_hypotesting_distancemetrics_oct2024.RDS')
	fm <- readRDS('resource_chp3/hypotesting_dredge_results/fishesmetric_glm_hypotesting_distancemetrics_MaxN_oct2024.RDS')
	pm <- readRDS('resource_chp3/hypotesting_dredge_results/largesharks_relativerisk_glm_hypotesting_distancemetricsoct2024.RDS')

# make table
	models <- c(sgm, fm, pm)
	
	sgf <- flextable::as_flextable(sgm)
	pf <- flextable::as_flextable(pm)
	ff <- flextable::as_flextable(fm)

	s1 <- tbl_regression(sgm,exp=FALSE,conf_level=0.95,label=list(st_d2sh='Dist. to Shore',st_dstc='Dist. to Central Mangroves',st_d2jetty='Dist. to Jetty','st_d2sh*st_dstc'='Dist. to Shore * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)
	f1 <- tbl_regression(fm, exp=FALSE,conf_level=0.95,label=list(st_d2sh='Dist. to Shore',st_d2jetty='Dist. to Nearest Jetty',st_mds='Mid-density seagrass',
		st_hds = 'Highest den. seagrass','st_d2jetty:st_d2sh'='Dist. to Nearest Jetty * Dist. to Shore'))%>%
		bold_p(t=0.05)	
	p1 <- tbl_regression(pm, exp=FALSE,conf_level=0.95, label=list(st_depth='Depth',st_d2shore='Dist. to Shore',st_distcmg='Dist. to Central Mangroves', 'st_depth:st_distcmg' = ' Depth * Dist. to Central Mangroves'))%>%
		bold_p(t=0.05)

	# side by side
	tbl_merge(tbls = list(s1,f1,p1), tab_spanner = c('Hypothesis 1', 'Hypothesis 2', 'Hypothesis 3'))
	# stacked 
	stacked <- tbl_stack(list(s1,f1,p1),group_header=c('1','2','3'))%>%
	bold_levels()
	
	show_header_names(stacked)
 
	responses <- as_tibble(x=c('Seagrass','Fish MaxN','Large Sharks'))%>%rename(Respones=value)

	stacked.summary<-stacked%>%
		modify_header(groupname_col = '**Hypothesis**',label='**Predictor**')%>%
		as_flex_table()%>%	
		autofit()

	## AIC of three. mdoels relative to each other
	AIC(sgm, fm, pm)

	save_as_image(stacked.summary,'resource_chp3/hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_Oct2024.png' ,webshot='webshot')
	save_as_docx(stacked.summary,path = 'resource_chp3/hypotesting_dredge_results/stackedsummarytables_hypotheses_chp3_seagrasses_fish_largesharks_glms_Oct2024.docx')

#############################
## - Residuals checking: histograms of residuals, and standardsised residuals versus fitted, looking for heteroscedasticity
#############################
    { #run to generate individual and panelled plots of best dredged GLMs of explanatory variables SG, fish and predation risk

    # calculate fitted and pull residuals 
		R1_sgm <- as_tibble(resid(sgm))%>%rename(resids = value)
		F1_sgm <- as_tibble(predict(sgm))%>%rename(fitted = value)
		diag_sgm <- bind_cols(F1_sgm, R1_sgm)

		R1_fm <- as_tibble(resid(fm))%>%rename(resids = value)
		F1_fm <- as_tibble(predict(fm))%>%rename(fitted = value)
		diag_fm <- bind_cols(F1_fm, R1_fm)

		R1_pm <- as_tibble(resid(pm))%>%rename(resids = value)
		F1_pm <- as_tibble(predict(pm))%>%rename(fitted = value)
		diag_pm <- bind_cols(F1_pm, R1_pm)

		# plots 
		resVfit_sgm <- ggplot(data = diag_sgm,aes(x= fitted, y = resids))+
			geom_point(col = '#3A6C74', fill = '#3A6C74')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Seagrasses')
		res_sgm_hist <- ggplot(data = R1_sgm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_sgm)/10000,col = '#3A6C74', fill = '#3A6C74')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_sgm_qq <- ggplot(data=F1_sgm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#3A6C74')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		resVfit_fm <- ggplot(data = diag_fm,aes(x= fitted, y = resids))+
			geom_point(col = '#2b0052', fill = '#2b0052')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Prey MaxN')
		res_fm_hist <- ggplot(data = R1_fm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_fm)/10000,col = '#2b0052', fill = '#2b0052')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_fm_qq <- ggplot(data=F1_fm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#2b0052')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		resVfit_pm <- ggplot(data = diag_pm,aes(x= fitted, y = resids))+
			geom_point(col = '#b27409', fill = '#b27409')+ 
			labs(subtitle = 'Residuals vs. Fitted', y = 'Residuals', x = 'Fitted')+
			theme_bw()+
			ggtitle('Predation Risk')
		res_pm_hist <- ggplot(data = R1_pm, aes(x = resids))+ 
			geom_histogram(binwidth = nrow(R1_pm)/1000,col = '#b27409', fill = '#b27409')+ 
			theme_bw()+
			labs(subtitle = 'Residuals', y = 'Frequency', x = 'Residuals')
		res_pm_qq <- ggplot(data=F1_pm, aes(sample = fitted))+
			stat_qq(size=2,pch=21,col = '#b27409')+
			labs(subtitle = 'Quantile-Quantile plot',x='Theoretical Quantiles',y='Standardised Residuals')+
			stat_qq_line(linetype=2, col='red')+
			theme_bw()

		diagnostics_sgm <- resVfit_sgm+res_sgm_hist+res_sgm_qq
		diagnostics_fm <- resVfit_fm+res_fm_hist+res_fm_qq
		diagnostics_pm <- resVfit_pm+res_pm_hist+res_pm_qq

		diagnostics_stacked <- diagnostics_sgm / diagnostics_fm / diagnostics_pm

		ggsave(diagnostics_sgm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_sagrasses_hypotestingGLM_oct2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)
		ggsave(diagnostics_fm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_fishMaxN_hypotestingGLM_oct2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)
		ggsave(diagnostics_pm, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_relrisk_hypotestingGLM_oct2024.png', device = 'png', unit = 'in', width = 8, height = 4, dpi = 850)		

		ggsave(diagnostics_stacked, file = 'resource_chp3/hypotesting_dredge_results/diagnostic_plots_hypotestingGLM_oct2024.png', device = 'png', unit = 'in', width = 8, height = 10, dpi = 850)

    }
