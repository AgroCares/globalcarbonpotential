# function to calculate the uncertainty on estimated C sequestration potential
# based on a Monte Carlo approach
# function is quite similar to fun_estimate_globa_cseq_potential.R

# require packages
require(sf);require(raster);require(data.table)

# clean environment
rm(list=ls())

# load data and ma-models
d1 <- fread('data/210322DB001.csv')
d1.ccID <-  fread('data/210322DB001_ccID.csv')
m1 <- as.data.table(readRDS('data/200609 model overview.rds'))

# save output in a list
out.list <- list()

# number of simulations
nsim = 1000

# set progressbar
pb <- txtProgressBar(title = "progress bar", max = nsim)

# do nsim simulations to estimate the global C sequestration potential
for (i in 1:nsim){
  
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
  
  # estimate effect of measure 1
  d2 <- merge(d1,m2[measure=='M1',.(wmean,wsd,climate)],by.x = 'climBeckGroup',by.y = 'climate')
  d2[,wmean2 := rnorm(.N,mean = wmean, sd = wsd),by = .(climBeckGroup)]
  d2[,wmean2 := pmax(0,wmean2)]
  d2[cat_fert=='NIFNOF',dc_m1 := wmean2]
  d2[cat_fert=='MIFNOF',dc_m1 := 0.5 * wmean2]
  d2[cat_fert=='MIFHOF',dc_m1 := 0.5 * wmean2]
  d2[cat_fert=='HIFNOF',dc_m1 := 0]
  d2[cat_fert=='HIFHOF',dc_m1 := 0]
  d2[,dc_tm1 := wmean2]
  d2[,c('wmean','wmean2','wsd') := NULL]
  
  # estimate effect of measure 2
  m2.sim <- m2[measure %in% c('M1','M2'),.(category,climate,wmean,wsd)]
  m2.sim[,wmean := rnorm(.N,mean = wmean, sd = wsd),by=.(climate,category)]
  m2.sim[,wmean := pmax(0,wmean)]
  m2.dcast <- dcast(m2.sim[,.(category,climate,wmean,wsd)],climate~category,value.var='wmean')
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
  m2.sim <- m2[measure=='M3' & category=='CC',.(category,climate,wmean,wsd)]
  m2.sim[,wmean := rnorm(.N,mean = wmean, sd = wsd),by=.(climate,category)]
  m2.sim[,wmean := pmax(0,wmean)]
  d2 <- merge(d2,m2.sim[,.(wmean,climate)],by.x = 'climBeckGroup',by.y = 'climate',all.x=TRUE)
  d2[cat_ccr == 'NOCATCH', dc_m3 := wmean]
  d2[cat_ccr == 'CATCH', dc_m3 := 0]
  d2[,wmean := NULL]
  
  # estimate effect of measure 4, tillage
  m2.sim <- m2[measure %in% c('M4'),.(category,climate,wmean,wsd)]
  m2.sim[,wmean := rnorm(.N,mean = wmean, sd = wsd),by=.(climate,category)]
  m2.sim[,wmean := pmax(0,wmean)]
  m2.dcast <- dcast(m2.sim[,.(category,climate,wmean)],climate~category,value.var='wmean')
  m2.dcast[,TILL_M := `NT-IT`]
  m2.dcast[,TILL_H := `NT-HT`]
  m2.dcast[,TILL_L := 0]
  m2.dcast <- rbind(m2.dcast,m2.dcast[climate=='oth'][,climate := 'trop'])
  m2.melt <- melt(m2.dcast,id.vars = 'climate',value.name = 'dc_m4',variable.name = 'tillage')
  m2.melt <- m2.melt[tillage %in% unique(d2$tillage)]
  d2 <- merge(d2,m2.melt,by.x = c('climBeckGroup','tillage'),by.y=c('climate','tillage'),all.x = TRUE)
  
  # estimate effect of measure 5, crop residue (there is only one estimate)
  m2.sim <- m2[measure %in% c('M5'),.(wmean,wsd)]
  m2.sim <- unique(m2.sim)
  m2.sim[,wmean := rnorm(.N,mean = wmean, sd = wsd)]
  m2.sim[,wmean := pmax(0,wmean)]
  m5.wmean = m2.sim[,wmean]
  d2[cat_crb == 'BURN',dc_m5 := m5.wmean]
  d2[cat_crb == 'INC', dc_m5 := 0]
  
  # in Mt (Megaton) (area in ha * ton C /ha / yr)
  scols <- c(paste0('dc_m',1:5),'dc_m2_corr')
  d3 <- d2[,lapply(.SD,function(x) round(sum(area * x)*0.001*0.001,0)),.SDcols = scols,by=.(climBeckGroup)]
  d3 <- rbind(d3,data.table(climBeckGroup='all',d3[,lapply(.SD,sum),.SDcols = scols]))
  d3[,id := 1]
  d4 <- melt(d3,id.vars=c('id','climBeckGroup'))
  d4 <- dcast(d4,id~climBeckGroup + variable)
  
  # save output in a list
  out.list[[i]] <- copy(d4)
  
  # print progress bar
  setTxtProgressBar(pb,i)
  
}

# combine all output into one data.table
out <- rbindlist(out.list)

# reformat the output
d5 <- melt(out,id.vars = 'id',measure=patterns("dc_m1$","dc_m2$",'dc_m3$','dc_m4$','dc_m5$','dc_m2_corr$'))
setnames(d5,c('id','climBeckGroup','dc_m1',"dc_m2",'dc_m3','dc_m4','dc_m5','dc_m2_corr'))
d5[climBeckGroup==1,climBeckGroup :='all']
d5[climBeckGroup==2,climBeckGroup :='oth']
d5[climBeckGroup==3,climBeckGroup :='stme']
d5[climBeckGroup==4,climBeckGroup :='temp']
d5[climBeckGroup==5,climBeckGroup :='trop']

# retreive mean and quantile information
cols <- c('dc_m1',"dc_m2",'dc_m3','dc_m4','dc_m5','dc_m2_corr')
d6 <-d5[ ,list(m1_mean = mean(dc_m1),
               m2_mean = mean(dc_m2),
               m3_mean = mean(dc_m3),
               m4_mean = mean(dc_m4),
               m5_mean = mean(dc_m5),
               m2cor_mean = mean(dc_m2_corr),
               m1_q95 = quantile(dc_m1,0.99),
               m2_q95 = quantile(dc_m2,0.99),
               m3_q95 = quantile(dc_m3,0.99),
               m4_q95 = quantile(dc_m4,0.99),
               m5_q95 = quantile(dc_m5,0.99),
               m2cor_q95 = quantile(dc_m2_corr,0.99),
               m1_q05 = quantile(dc_m1,0.01),
               m2_q05 = quantile(dc_m2,0.01),
               m3_q05 = quantile(dc_m3,0.01),
               m4_q05 = quantile(dc_m4,0.01),
               m5_q05 = quantile(dc_m5,0.01),
               m2cor_q05 = quantile(dc_m2_corr,0.01)
),by = .(climBeckGroup)]

# reformat the output per measure
d6 <- melt(d6,id.vars='climBeckGroup',measure=patterns('^m1_','^m2_','^m3_','^m4_','^m5_','^m2cor_'))
d6[variable==1,variable := 'mean']
d6[variable==2,variable := 'q95']
d6[variable==3,variable := 'q05']
setorder(d6,'climBeckGroup')
setnames(d6,c('climBeckGroup','stats',c('m1','m2','m3','m4','m5','m2cor')))
cols <- c('m1','m2','m3','m4','m5','m2cor')
d6[,c(cols) := lapply(.SD,function(x) round(x,0)),.SDcols = cols]

d7 <- melt(d6,id.vars=c('climBeckGroup','stats'))
d7 <- dcast(d7,climBeckGroup+variable~stats)
d7[,tab := paste0(mean," (",q05," - ",q95,")")]
d7 <- dcast(d7,climBeckGroup~variable,value.var = 'tab')
d7[c(3,4,5,2,1)]
