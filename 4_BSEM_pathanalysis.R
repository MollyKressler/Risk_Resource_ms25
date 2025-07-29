## Bayesian Structural Equartion Model & Path analysis for MS 'Bayesian Structural Equation Models reveal risk and resource trade-offs in a juvenile elasmobranch'

## Code and approach here heavily inspired by modelling code from paper Bennett, S., Harris, M. P., Wanless, S., Green, J. A., Newell, M. A., Searle, K. R., & Daunt, F. (2022). Earlier and more frequent occupation of breeding sites during the non‐breeding season increases breeding success in a colonial seabird. Ecology and Evolution, 12(9). https://doi.org/10.1002/ece3.9213. Thank you to S Bennett for sharing their JAGs code with Rich B Sherley. 

## created by Molly M Kressler 


###################################
########## RUN AT OPEN ############
###################################

## Load Workspace, local macbook

pacman::p_load(MCMCvis,tidyverse,sf,nimble,devtools,flextable,arm,webshot2,sfdep,sp,spdep,beepr,HDInterval, patchwork, cowplot, tidybayes)
setwd('/Users/mollykressler/Documents/Documents - Molly’s MacBook Pro/data_phd')


###################################
##########     END     ############
###################################

#####################################################
###### DF 1, for process model for sharkiness  ######
### Define Sharkiness
## counts of detections per individual per Site (buffer/receiver)

pointdata<-read.csv('pointdata_juvlemons_withAllCoV_MKthesis20192020.csv')%>%
    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
    mutate(buffIDnum = parse_number(buffID))
stopifnot(nrow(pointdata)==560*2) # check 
sapply(pointdata,class)
summary(pointdata)
  
 #add relative proportion of detections of juveniles to path inference figure.
    #summarise across tidal phases (L/H)
    names(pointdata)
    point.sf <- st_as_sf(st_read('pointdata_juvlemons_withAllCoV_MKthesis20192020.shp'),crs = 'WGS84')
    total = sum(point.sf$n_juv)
    sum.dett <- point.sf %>% 
      dplyr::select(buffID, n_juv, tidephs,geometry) %>%
      group_by(buffID)%>%
      summarise(n = sum(n_juv))%>%
      st_as_sf()
    sum.dett

####################################################
###### DF 2, for process model for fishiness  ######

hexdata<-read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%
    mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))
stopifnot(nrow(hexdata)==5296) # check 
summary(hexdata)

hexsf <- st_as_sf(st_read('hexdata_juvLemonsPrUse_withAllCov_MKthesis.shp'),crs='WGS84')%>%
  mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
  mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
  rename(standard.hexshark = st_shark,
    standard.hexfish = st_maxN,
    standard.hexdist2shore = st_d2sh,
    standard.hexdistcmg = st_dstc,
    standard.hexmds = st_mds,
    standard.hexhds = st_hds,
    standard.hexdist2jetty = st_d2jetty)


##################################################################
##################################################################

##################### 
## - MODEL6

  modelCode6<-nimbleCode({
      

    ##########################
    ######### priors #########
    
    #### prior for sharkiness #### 
    for(i in 1:7){
       a[i] ~ dnorm(0,.01) 
      }
    for(i in 1:2){
       g[i] ~ dnorm(0,.01)# for pred and fish only models, for interpretting coefficients 
      }
      # prior for residual variance - sharks 
      tau.shark ~ dgamma(0.01,0.01) 
      sigma.shark <- sqrt(1/tau.shark)
      tau.shark.hex ~ dgamma(0.01, 0.01) # for glm of fish on sharks, at hexagon level
      sigma.shark.hex <- sqrt(1/tau.shark.hex)
    
    #### prior for predators #### 
    for(i in 1:5){
       e[i] ~ dnorm(0,.001) 
      }
      # prior for intercept - sharks 
      f ~ dnorm(0,.001)
      # prior for residual variance - sharks 
      tau.pred ~ dgamma(0.01,0.01) # 
      sigma.pred <- sqrt(1/tau.pred)
    
    # group-level effect of buffID on sharkiness and predators.
    for(i in 1:B){
      epsi_shark[i]~ dnorm(0,tau.epsi_shark)
    }
      # prior for intercept - sharks 
      b ~ dnorm(0,.001)
      tau.epsi_shark ~ dgamma(0.01,0.01)
      sigma.epsi_shark <- sqrt(1/tau.epsi_shark)
    
    #### prior for fishiness (cross metric, hexagons) ####
    for(i in 1:6){
        c[i] ~ dnorm(0,0.001) 
      }
      # prior for intercept - fishiness
      d ~ dnorm(0,0.001)
      
      # prior for residual variance (precision) - fishiness
      tau.fish ~ dgamma(0.01,0.01) 
      sigma.fish <- sqrt(1/tau.fish) # gives sd 
    
    #### prior for seagrasses PCA @ hexagon level ####
    for(i in 1:4){
      j[i] ~ dnorm(0,0.01) 
          }
      # prior for intercept - hex 
      k ~ dnorm(0,0.001) # low
          
      # prior for residual variance (precision) - hex sg PCA 
      tau.hexsg ~ dgamma(0.01,0.01)
      sigma.hexsg <- sqrt(1/tau.hexsg)
  

    #########################################################
    ######### Likelihoods - data and process models ######### 
     ## informed by hypothesis exploration 

    ### data model for seagrasses - hexagons
      for(i in 1:hex.N){
      standard.hexsg[i] ~ dnorm(hexsg.mu[i],tau.hexsg)

      hexsg.mu[i] <- k + j[1]*standard.hexdist2shore[i] + j[2]*standard.hexdistcmg[i] + j[3]*standard.hexdist2jetty[i] + j[4]*standard.hexdist2shore[i]*standard.hexdistcmg[i] 
    }
    
    ### data model for fishiness - hexagons 
     for(i in 1:hex.N){
      standard.hexfish[i] ~ dnorm(fish.mu[i],tau.fish) 

      fish.mu[i] <- d + c[1]*standard.hexdist2jetty[i] + c[2]*standard.hexdist2shore[i] + c[3]*standard.hexmds[i]+ c[4]*standard.hexhds[i]+ c[5]*standard.hexdist2jetty[i]*standard.hexdist2shore[i] + c[6]*hextide[i]

      # glm for the effect of fish (only)
      standard.hexshark[i] ~ dnorm(fishonly.mu[i], tau.shark.hex)
      fishonly.mu[i] <- d + g[1]*standard.hexfish[i] 
    }
  
    ### data model for large shark detectons - point data 
     for(i in 1:point.N){
      zlogit.sqzrisk[i] ~ dnorm(pred.mu[i], tau.pred)

      pred.mu[i] <- f + e[1]*standard.pointdist2shore[i] + e[2]*standard.pointdepth[i] + e[3]*standard.pointdistcmg[i] + e[4]*standard.pointdepth[i]*standard.pointdistcmg[i]+ e[5]*pointtide[i]


      ## glm for the effect of predation (only) 
      standard.shark3[i] ~ dnorm(predonly.mu[i], tau.shark)
      predonly.mu[i] <-  epsi_shark[buffID[i]] + g[2]*zlogit.sqzrisk.pred[i]
     }

    ### data model for sharkiness - pointdata
     for(i in 1:point.N){
      standard.shark[i] ~ dnorm(shark.mu[i],tau.shark)   
      shark.mu[i] <- epsi_shark[buffID[i]] + a[1]*standard.fish.pred[i] + a[2]*standard.pointdist2shore[i] + a[3]*standard.pointdistcmg[i] + a[4]*standard.sg[i] + a[5]*standard.pointdepth[i]+ a[6]*standard.pointdist2jetty[i]+ a[7]*zlogit.sqzrisk[i] 

      }  

    ######### Derived Parameters #########
    # for estmating total pathways 
    ## coefficients for distance metrics are from process models of those predictors. leave out coefficients for interactions if the fixed effects coefficients are already included.

    path[1] <-  a[4] * j[1] * j[2] * j[3] # seagrass and abiotics

    path[2] <-  a[1] * c[1] * c[2] * c[3] * c[4] * c[6] # fish, the things that effect fish, and tide

    path[3] <-  a[7] * e[1] * e[2] * e[3] * e[5] # predator risk, the things that effect risk, and tide

    ## need to calculate the values from the abiotics along the paths to the initial parameter, for ease

    value[1] <-  c[1] * c[2] * c[3] * c[4]   # path 2, abiotics
    value[2] <- e[1] * e[2] * e[3]  # path 3, abiotics

  }) # end of model code 


    ##########################
    ## Compile the model code
    ##########################
    
    myConstants6<-list(point.N=1120,hex.N=5296,B=max(pointdata$buffIDnum),buffID=pointdata$buffIDnum)
    
    myData6<-list(
      # tell nimble the covariates 
      # point data
      standard.shark = pointdata$st_shark,
      standard.pointdist2shore = pointdata$st_d2shore,
      standard.pointdistcmg = pointdata$st_distcmg,
      standard.fish.pred = pointdata$st_maxN,
      standard.sg = pointdata$st_SG,
      standard.pointdepth= pointdata$st_depth,
      standard.pointdist2jetty= pointdata$st_d2jetty,
      zlogit.sqzrisk = pointdata$st_risk,
      zlogit.sqzrisk.pred = pointdata$st_risk,
      standard.shark3 = pointdata$st_shark,
      pointtide = pointdata$tide,
      # hex data 
      standard.hexfish = hexdata$st_maxN,
      standard.hexdist2shore = hexdata$st_d2sh,
      standard.hexdistcmg = hexdata$st_dstc,
      standard.hexdist2jetty = hexdata$st_d2jetty,
      standard.hexmds = hexdata$st_mds,
      standard.hexhds = hexdata$st_hds,
      standard.hexshark = hexdata$st_shark,
      hextide = hexdata$tide,
      standard.hexsg = hexdata$st_SG
      )
    
    init.values6<-list(a=rnorm(7,0,1),
                        b=rnorm(1),
                        c=rnorm(6,0,1),
                        d=rnorm(1),
                        e=rnorm(5,0,1),
                        f=rnorm(1),
                        j=rnorm(4,0,1),
                        k=rnorm(1),
                        g=rnorm(2,0,1),
                        path=rnorm(3,0,1),
                        value=rnorm(2,0,.05),
                        epsi_shark=rgamma(35,1,0.1),
                        tau.shark=rgamma(1,0.01,0.01),
                        tau.fish=rgamma(1,0.01,0.01),
                        tau.pred=rgamma(1,0.01,0.01),
                        tau.epsi_shark=rgamma(1,0.01,0.01),
                        tau.hexsg=rgamma(1,0.01,0.01),
                        tau.shark.hex = rgamma(1, 0.01, 0.01)
    )
    
    ## model4 define and compile
    model6<-nimbleModel(code=modelCode6, name="model6",data=myData6,constants = myConstants6,inits=init.values6) #define the model
    
    model6$calculate() # if = NA, indicates missing or invalid initial values, and you have to fix the model until it is numeric.
    model6$initializeInfo()
    
    Cm6<-compileNimble(model6) # compile the model

    
    ##########################
    ## Compile & Run the MCMC
    ##########################
    
    ## model5
    conf6<- configureMCMC(model6,monitors=c('a','b','c','d','j','k','e','f','g','tau.epsi_shark','tau.fish','tau.shark','tau.pred','tau.hexsg','path','value'),onlySlice=FALSE) 
    MCMC_model6 <- buildMCMC(conf6,na.rm=TRUE)
    ccMCMC6 <-compileNimble(MCMC_model6, project = model6)
    samples6 <- runMCMC(ccMCMC6,niter=10000, nburnin=2000, nchains=3,samplesAsCodaMCMC = TRUE) 

    summary(samples6)
    
    saveRDS(samples6,'resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.RDS')

    # trace and density plots
    
    MCMCtrace(samples6,pdf=TRUE,ind=TRUE, Rhat=TRUE, n.eff=TRUE) # ind = TRUE, separate density lines per chain. # pdf = FALSE, don't export to a pdf automatically. 
    
     ###########################################################
    ## Summary Table & Caterpillar plots with MCMCvis & tidybayes to show small values ##
    ###########################################################
    
    # import RDS, local macbook 
    samplesList6a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.RDS')

    mcmc_summary_Cmodel6_samplesListfromRDS<-MCMCsummary(samplesList6a,round=3,probs=c(0.05,0.95),pg0=TRUE)%>%
      tibble::rownames_to_column()%>%
      rename_with(str_to_title)%>%
      rename('pg0'='P>0')%>%
      mutate(pg00 = case_when(Mean >= 0 ~ as.numeric(pg0), Mean < 0 ~ 1-as.numeric(pg0), .default = as.numeric(pg0)))%>%
      rename(Parameter = Rowname, 'Prop. of posterior with \n\ same sign as estimate' = 'pg00', Estimate = 'Mean','lower'='5%',upper='95%')%>%
      mutate(CI = paste0('[',lower,',',upper,']'),.after='Estimate')%>%
      dplyr::select(-lower,-upper,-Sd, -pg0)%>% 
      filter(Parameter!='value[1]', Parameter!='value[2]',Parameter!='value[3]')%>%
      flextable()%>%
      theme_zebra()%>%
      set_header_labels(rowname = 'Coefficient',SD='Sd')%>%
      align(align = 'center', part = 'all')%>%
      font(fontname = 'Arial', part = 'all')%>%
      color(color='black',part='all')%>%
      fontsize(size = 10, part = 'all')%>%
      autofit()
    mcmc_summary_Cmodel6_samplesListfromRDS
    
    save_as_image(mcmc_summary_Cmodel6_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples__model6_niter5000_burn1000_chains3_jan2025.png',res=850)  
    save_as_docx(mcmc_summary_Cmodel6_samplesListfromRDS,path='resource_chp3/nimblemodel_outputs/mcmcsamples__model6_niter5000_burn1000_chains3_jan2025.docx')  
    ## for chp3 ms v6 

      save_as_image(mcmc_summary_Cmodel6_samplesListfromRDS,path='chp3_ms_v6_figures_tables/mcmcsamples__model6_niter5000_burn1000_chains3_jan2025.png',res=850)  
    save_as_docx(mcmc_summary_Cmodel6_samplesListfromRDS,path='chp3_ms_v6_figures_tables/mcmcsamples__model6_niter5000_burn1000_chains3_jan2025.docx')  

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
    

##################################################
## Inference from Paths ##
##################################################

    ## For local macbook
    samplesList4 <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model4_niter20000_burn12000_chains3_4may2024.RDS')
    
    samplesList5a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model5a_niter5000_burn1000_chains3_July2024.RDS')

    samplesList6a <- readRDS('resource_chp3/nimblemodel_outputs/mcmcsamples_model6_niter5000_burn1000_chains3_jan2025.RDS')

    pointdata<-read.csv('pointdata_juvlemons_withAllCoV_MKthesis20192020.csv')%>%
    mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
    mutate(buffIDnum = parse_number(buffID))

    hexdata <- read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%
      mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
      mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
      rename(
        standard.hexdist2shore = st_d2sh, 
        standard.hexdepth = st_depth,
        standard.hexdistcmg = st_dstc,
        hextide = tide,
        zlogit.sqzrisk = st_risk
        )
    hexsf <- st_as_sf(st_read('hexdata_juvLemonsPrUse_withAllCov_MKthesis.shp'),crs='WGS84')%>%
      mutate(st_shark = (m5_PUSE - mean(m5_PUSE))/sd(m5_PUSE))%>%
      mutate(tide = case_when(tidephs == 'L' ~ 1, tidephs == 'H' ~ 0))%>%
      rename(
        standard.hexdist2shore = st_d2sh, 
        standard.hexdepth = st_depth,
        standard.hexdistcmg = st_dstc,
        hextide = tide,
        zlogit.sqzrisk = st_risk
        )
      nrow(hexsf) # should be 5296
    
  ## Make test sample data frames - based on the hexagon df

    hexsamp <- hexdata %>%
      sample_n(5)%>%
      as.data.frame() 


  ## From model4 samplesList, make objects with all draws of each coeefficient from path 2 and path 3, e.g. j4. 
    # path 3: e1, e2, e3, e4, e5, a7
    e1.ch <-c(samplesList6a$chain1[,16], samplesList6a$chain2[,16],samplesList6a$chain3[,16] )
    e2.ch <-c(samplesList6a$chain1[,17], samplesList6a$chain2[,17],samplesList6a$chain3[,17] )
    e3.ch <-c(samplesList6a$chain1[,18], samplesList6a$chain2[,18],samplesList6a$chain3[,18] )
    e5.ch <-c(samplesList6a$chain1[,20], samplesList6a$chain2[,20],samplesList6a$chain3[,20] )
    a7.ch <-c(samplesList6a$chain1[,7], samplesList6a$chain2[,7],samplesList6a$chain3[,7] )
  ## Run loops for each path, to calculate path estimate given a hexagon cell

      ## Define size of objects/matrices
      J = 3 * 8000 # chains x iterations, 8000 in June, 4000 in July, 8000 in january
      n.hex = nrow(hexdata) 

      ## set up prediction df

        preds.path3 = as.data.frame(matrix(NA,ncol = J, nrow = n.hex))
        head(preds.path3)


      ## write progress bar function
        pb <- txtProgressBar(min = 1, max = n.hex, style = 3)

      ## path 3 loop
      for(i in 1:n.hex){
        setTxtProgressBar(pb, i)
        for(j in 1:J){
          preds.path3[i,j] <-  (e1.ch[j]*hexdata$standard.hexdist2shore[i]) + (e2.ch[j]*hexdata$standard.hexdepth[i]) + (e3.ch[j]*hexdata$standard.hexdistcmg[i]) + (e5.ch[j]*hexdata$hextide[i]) + (a7.ch[j]*hexdata$zlogit.sqzrisk[i]) 
        }
      }

      ## save calculations
        saveRDS(preds.path3,'resource_chp3/path_inference/path3_estimates_at_hexagons_model6aJuly2024_calcJan2025.RData')

  ## Calculate marginal means, and HDI (highest density intervals)

    ## deprecated and inaccurate way to calculate HDIs - HDInterval::hdi needs values in columns. 

      mean(as.numeric(preds.path3[1,]))

      cols <- c('jcode','mean', 'lower', 'upper')
      p2pred <- as.data.frame(matrix(ncol=4, nrow = n.hex))
      p3pred <- as.data.frame(matrix(ncol=4, nrow = n.hex))
      colnames(p2pred) = cols
      colnames(p3pred) = cols
      p2pred$jcode <- as.character(hexdata$jcode)
      p3pred$jcode <- as.character(hexdata$jcode)
      head(p2pred)

      pb <- txtProgressBar(min = 1, max = n.hex, style = 3)
      
      for(i in 1:n.hex){
        p2pred[i,2] <- mean(as.numeric(preds.path2[i,]))
        p2pred[i,3] <- hdi(preds.path2[i,])[2]
        p2pred[i,4] <- hdi(preds.path2[i,])[1]
        setTxtProgressBar(pb, i)
        p3pred[i,2] <- mean(as.numeric(preds.path3[i,]))
        p3pred[i,3] <- hdi(preds.path3[i,])[2]
        p3pred[i,4] <- hdi(preds.path3[i,])[1]
      };beep(3)

      head(p3pred)

    ## Updated approach: 5 March 2024
    
    p3 <- readRDS('resource_chp3/path_inference/path3_estimates_at_hexagons_model6aJuly2024_calcJan2025.RData')
    
    hexdata <- read.csv('hexdata_juvLemonsPrUse_withAllCov_MKthesis.csv')%>%mutate(jcode=as.character(jcode))
      head(hexdata)

    ## re-format data - but tides! predictiosn follow the hexdata index order - which puts all LOW tide for each hexagon, then all high tide. So it will be 2648 rows of LOW + jcode, then 2648 rows of HIGH + jcode
    dat.p3 <- as.matrix(p3)
    out.p3 <- data.frame(
        Index = rep(0,dim(dat.p3)[1]),
        Mean = rep(0,dim(dat.p3)[1]),
        Low = rep(0,dim(dat.p3)[1]),
        Upp = rep(0,dim(dat.p3)[1]))
    head(out.p3)

    ## calculate means and HDIs for each hexagon at each tidephs
    
    for(i in 1:dim(dat.p3)[1]){
        out.p3$Index[i] = i
        out.p3$Mean[i] = mean(dat.p3[i,])
        out.p3$Low[i] = HDInterval::hdi(dat.p3[i,])[1]
        out.p3$Upp[i] = HDInterval::hdi(dat.p3[i,])[2]
        } # takes a beat
      head(out.p3)

      ## attach jcodes back 
      out.p3 <- out.p3 %>%
          mutate(jcode = hexdata$jcode, tidephs = hexdata$tidephs, .before = 'Index')%>%
        dplyr::select(-Index)
      head(out.p3)

  ## Save path estimates and path means + HDI dfs

    saveRDS(out.p3, 'resource_chp3/path_inference/path3_means_andHDI_at_hexagons_model6_jan2025.RData')







