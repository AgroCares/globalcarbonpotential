# script to merge meta-analytical model to global dataset to get an estimate of the potential C sequestration
# source data has been prepared in ppr_datasets_for_upscaling.R
# data has been merged into one global set with all relevant properties for upscaling in ppr_datasets_merging.R
# this script gives the merging of the global set with the meta-analytical estimates derived from the meta-analytical models collected.

# require packages
require(sf);require(raster);require(data.table)

# load in the data --------

  # clean environment
  rm(list=ls())
  
  # load the database
  d1 <- fread('data/210322DB001.csv')
  d1.ccID <-  fread('data/210322DB001_ccID.csv')
  
  # load the meta-analytcial models
  m1 <- as.data.table(readRDS('data/200609 model overview.rds'))

# estimate global C sequestration potential ---------

  # calculate weighted mean effect per measure for categories and climate zones
  m1.full <- m1[,.(mean = round(mean(Ct_yr),3),
                   wmean = round(weighted.mean(Ct_yr,Ct_ivar),3),
                   wsd = 1/sqrt(sum(1/Ct_yr_se^2))),by=.(measure,category,climate)]
  
  # calculate weighted mean effect per measure for categories without climate zones
  m1[,Ct_ivar_global := 1/((Ct_yr_se^2)*n),by=.(measure,category)]
  m1.global <- m1[,.(mean = round(mean(Ct_yr),3),
                     wmean = round(weighted.mean(Ct_yr,Ct_ivar_global),3),
                     wsd = 1/sqrt(sum(1/Ct_yr_se^2))),by=.(measure,category)]
  
  # combine both model estimates
  m1.global[,climate := 'oth']
  m2 <- rbind(m1.full,m1.global)
  
  # reformat table for reporting purposes (NOT USED)
  m3 <- copy(m2)
  m3[,value := paste0(round(wmean,2)," ± ",round(wsd,2))]
  d3 <- dcast(m3,measure + category~ climate, value.var = 'value',fill=' ')
  setcolorder(d3,c('measure','category','temp','oth'))
  
  
  # estimate effect of measure 1
  d2 <- merge(d1,m2[measure=='M1',.(wmean,climate)],by.x = 'climBeckGroup',by.y = 'climate')
  d2[cat_fert=='NIFNOF',dc_m1 := wmean]
  d2[cat_fert=='MIFNOF',dc_m1 := 0.5 * wmean]
  d2[cat_fert=='MIFHOF',dc_m1 := 0.5 * wmean]
  d2[cat_fert=='HIFNOF',dc_m1 := 0]
  d2[cat_fert=='HIFHOF',dc_m1 := 0]
  d2[,dc_tm1 := wmean]
  d2[,wmean := NULL]
  
  # estimate effect of measure 2
  m2.dcast <- dcast(m2[measure %in% c('M1','M2'),.(category,climate,wmean)],climate~category,value.var='wmean')
  m2.dcast[,NIFNOF := `OF-NF`]
  m2.dcast[,MIFNOF := `COF-NF` - `IF-NF`]
  m2.dcast[,MIFHOF := `COF-NF` - `OF-NF`]
  m2.dcast[,HIFNOF := 0.5 * MIFNOF]
  m2.dcast[,HIFHOF := 0]
  m2.melt <- melt(m2.dcast,id.vars = 'climate',value.name = 'dc_m2',variable.name = 'cat_fert')
  d2 <- merge(d2,m2.melt,by.x = c('climBeckGroup','cat_fert'),by.y=c('climate','cat_fert'),all.x = TRUE)
  
  # do correction for amount of manure available (CN or 10, humification coefficient of 50%)
  d2[,dc_m2_tonC := nloss_kgha_mean * 10 * 0.001 * 0.5]
  d2[,dc_m2_corr := pmin(dc_m2,dc_m2_tonC)]
  
  # estimate effect of measure 3 cover crops (3a)
  d2 <- merge(d2,m2[measure=='M3' & category=='CC',.(wmean,climate)],by.x = 'climBeckGroup',by.y = 'climate',all.x=TRUE)
  d2[cat_ccr == 'NOCATCH', dc_m3 := wmean]
  d2[cat_ccr == 'CATCH', dc_m3 := 0]
  d2[,wmean := NULL]
  
  # estimate effect of measure 4, tillage
  m2.dcast <- dcast(m2[measure %in% c('M4'),.(category,climate,wmean)],climate~category,value.var='wmean')
  m2.dcast[,TILL_M := `NT-IT`]
  m2.dcast[,TILL_H := `NT-HT`]
  m2.dcast[,TILL_L := 0]
  m2.dcast <- rbind(m2.dcast,m2.dcast[climate=='oth'][,climate := 'trop'])
  m2.melt <- melt(m2.dcast,id.vars = 'climate',value.name = 'dc_m4',variable.name = 'tillage')
  m2.melt <- m2.melt[tillage %in% unique(d2$tillage)]
  d2 <- merge(d2,m2.melt,by.x = c('climBeckGroup','tillage'),by.y=c('climate','tillage'),all.x = TRUE)
  
  # estimate effect of measure 5, crop residue (there is only one estimate)
  m5.wmean = unique(m2[measure=='M5',.(wmean)])
  d2[cat_crb == 'BURN',dc_m5 := m5.wmean]
  d2[cat_crb == 'INC', dc_m5 := 0]
  
  # get totals for reporting analysis for the paper
  d2[,round(sum(area) *100/sum(d2$area)),by=.(climBeckGroup)]
  d2[,round(sum(area) *100/sum(d2$area)),by=.(tillage)]
  d2[,round(sum(area) *100/sum(d2$area)),by=.(cat_fert)]
  d2[,round(sum(area) *100/sum(d2$area)),by=.(cat_crb)]
  d2[,round(sum(area) *100/sum(d2$area)),by=.(cat_ccr)]
  
  # estimate impact in Mt (Megaton) (area in ha * ton C /ha / yr)
  scols <- c(paste0('dc_m',1:5),'dc_m2_corr')
  d3 <- d2[,lapply(.SD,function(x) round(sum(area * x)*0.001*0.001,0)),.SDcols = scols,by=.(climBeckGroup)]
  d3 <- rbind(d3,data.table(climBeckGroup='all',d3[,lapply(.SD,sum),.SDcols = scols]))
  d3[c(2:4,1,5)]


# - make spatial file for mapping ---------------

  # estimate effect of measure 1
  d2 <- merge(d1.ccID,m2[measure=='M1',.(wmean,climate)],by.x = 'climBeckGroup',by.y = 'climate')
  d2[cat_fert=='NIFNOF',dc_m1 := wmean]
  d2[cat_fert=='MIFNOF',dc_m1 := 0.5 * wmean]
  d2[cat_fert=='MIFHOF',dc_m1 := 0.5 * wmean]
  d2[cat_fert=='HIFNOF',dc_m1 := 0]
  d2[cat_fert=='HIFHOF',dc_m1 := 0]
  d2[,wmean := NULL]
  
  # estimate effect of measure 2
  m2.dcast <- dcast(m2[measure %in% c('M1','M2'),.(category,climate,wmean)],climate~category,value.var='wmean')
  m2.dcast[,NIFNOF := `OF-NF`]
  m2.dcast[,MIFNOF := `COF-NF` - `IF-NF`]
  m2.dcast[,MIFHOF := `COF-NF` - `OF-NF`]
  m2.dcast[,HIFNOF := 0.5 * MIFNOF]
  m2.dcast[,HIFHOF := 0]
  m2.melt <- melt(m2.dcast,id.vars = 'climate',value.name = 'dc_m2',variable.name = 'cat_fert')
  d2 <- merge(d2,m2.melt,by.x = c('climBeckGroup','cat_fert'),by.y=c('climate','cat_fert'),all.x = TRUE)
  
  # do correction for amount of manure available (CN or 10, humification coefficient of 50%)
  d2[,dc_m2_tonC := nloss_kgha * 10 * 0.001 * 0.5]
  d2[,dc_m2_corr := pmin(dc_m2,dc_m2_tonC)]
  
  # estimate effect of measure 3 cover crops (3a)
  d2 <- merge(d2,m2[measure=='M3' & category=='CC',.(wmean,climate)],by.x = 'climBeckGroup',by.y = 'climate',all.x=TRUE)
  d2[cat_ccr == 'NOCATCH', dc_m3 := wmean]
  d2[cat_ccr == 'CATCH', dc_m3 := 0]
  d2[,wmean := NULL]
  
  # estimate effect of measure 4, tillage
  m2.dcast <- dcast(m2[measure %in% c('M4'),.(category,climate,wmean)],climate~category,value.var='wmean')
  m2.dcast[,TILL_M := `NT-IT`]
  m2.dcast[,TILL_H := `NT-HT`]
  m2.dcast[,TILL_L := 0]
  m2.dcast <- rbind(m2.dcast,m2.dcast[climate=='oth'][,climate := 'trop'])
  m2.melt <- melt(m2.dcast,id.vars = 'climate',value.name = 'dc_m4',variable.name = 'tillage')
  m2.melt <- m2.melt[tillage %in% unique(d2$tillage)]
  d2 <- merge(d2,m2.melt,by.x = c('climBeckGroup','tillage'),by.y=c('climate','tillage'),all.x = TRUE)
  
  # estimate effect of measure 5, crop residue (there is only one estimate)
  m5.wmean = unique(m2[measure=='M5',.(wmean)])
  d2[cat_crb == 'BURN',dc_m5 := m5.wmean]
  d2[cat_crb == 'INC', dc_m5 := 0]
  
  # in Mt (Megaton) (area in ha * ton C /ha / yr)
  d3 <- d2[,lapply(.SD,function(x) round(sum(area_till_cr * x)*0.001*0.001,5)),.SDcols = c(paste0('dc_m',1:5),'dc_m2_corr'),by=.(ccID)]
  
  # in tonnes (area in ha * ton C /ha / yr)
  d3.a <- d2[,lapply(.SD,function(x) round(sum(area_till_cr * x),5)),.SDcols = c(paste0('dc_m',1:5),'dc_m2_corr'),by=.(ccID)]
  d3.b <- d2[,list(area = sum(area_till_cr)),by = .(ccID)]
  d4 <- merge(d3.a,d3.b,by='ccID')
  cols <-  c(paste0('dc_m',1:5),'dc_m2_corr')
  # convert to kg / ha, mean value per ccID
  d4[,c(cols) := lapply(.SD,function(x) round(x * 1000 / area,1)),.SDcols = cols]
  
  # read in shape per country-climate
  cc <- readRDS('../data/products/countryclimate.rds')
  cc$area_poly = st_area(cc)/10000
  cc <- as.data.table(cc)
  cc2 <- merge(cc,d4[,mget(c('ccID',c(paste0('dc_m',1:5),'dc_m2_corr')))],by.x = 'ID',by.y='ccID',all.x=TRUE)
  cc2[, area_poly := as.numeric(area_poly)]
  
  # update NA cells with median value per county / climate code (GRID-cell)
  cc2[,area_cor := as.numeric(area_poly / sum(area_poly[!is.na(dc_m1)])),by =.(name0,GRIDCODE)]
  cols = c(paste0('dc_m',1:5),'dc_m2_corr')
  cc2[,c(paste0(cols,'_median')) := lapply(.SD,function(x) fifelse(is.na(x),sum(x * area_cor,na.rm=T),x)),.SDcols = cols,by=.(name0,GRIDCODE)]
  
  # update NA cells with zero
  cc2[,c(cols) := lapply(.SD,function(x) fifelse(is.na(x),0,x)),.SDcols = cols]
  
  # convert to sf
  cc.sf <- st_as_sf(cc2)
  cc.sf <- st_cast(cc.sf,'MULTIPOLYGON')
  
  # save as gpkg file (needed for figure 1 in article)
  sf::st_write(cc.sf,'products/carbon_seq_pot_v4.gpkg')  
