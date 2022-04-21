aoi#install.packages('reader')
library(reader)
library(rgdal)
library(tidyverse)
library(zoo)
library(lubridate)
require(hydroGOF)
# swat folder
maindir = 'D:/sda'
basinNames = c('vye','gvo','aho','bye','slu', 'chu','gso','nkh','xla')
riverSystems = c('ht','sk', 'sk','mk','sk','ht','mk','ca','ma')

figPath = file.path(maindir,'0_figures')

#own-built function
source('C:/Users/Admin/Dropbox/RStudy/0.Code/built_in/sm_support_functions.r')
source('C:/Users/Admin/Dropbox/RStudy/0.Code/built_in/plot_support_functions.r')

dateAll = seq(as.Date('01-01-2015', format = '%m-%d-%Y'),
           as.Date('12-31-2019', format = '%m-%d-%Y'), by = 'day')
nb = length(basinNames)
nwu = 365 + 366 # number of day for the year of 2015 and 2016



for(ib in 3:3){
  folderName = paste('swat','_',riverSystems[ib],'_',
                     basinNames[ib], sep = '')
  basinShp = readOGR(paste(maindir,'/',folderName,'/','00dataPreparation','/','knownSubbasins','/',
                           basinNames[ib],'_subs.shp', sep = ''))
  nSubs = length(basinShp)
  
  # read and process the raw SMAP data for each HRU -------------------------
  rawSM = read.csv(file.path(maindir, folderName,'00dataPreparation','processedSM','smSMAP.csv'))
  nsim = nrow(rawSM)



  
  # read sub-basins and hrus information
  temp = read_fwf(file.path(maindir, folderName, '03openloop', 'Project','input.std'), skip = 39,
                  n_max = nSubs+1)
  hruInfo = temp[2:nrow(temp),1:4]
  colnames(hruInfo)  = temp[1,1:4]
  hruInfo = data.frame(sub = as.numeric(hruInfo$Sub),
                       hru = as.numeric(hruInfo$`#HRUs`),
                       lb = numeric(nrow(hruInfo)),
                       ub = numeric(nrow(hruInfo)))
  
  # indentify lower bound and upper bound
  hruInfo$lb[1] = hruInfo$sub[1]
  hruInfo$ub[1] = hruInfo$hru[1]
  
  for(ii in 2:nrow(hruInfo)){
    hruInfo$lb[ii] = hruInfo$ub[ii-1]+1
    hruInfo$ub[ii] = hruInfo$lb[ii] + hruInfo$hru[ii] -1
  }
  nhru =max(hruInfo$ub)
  # dertermine column order of rawSM
  rawOrder = colnames(rawSM)[-1]
  rawOrdernum = as.numeric(str_remove(rawOrder,'sm'))
  
  
  sm4da = mat.or.vec(nsim, nhru)
  for(isub in 1:nSubs){
    id = which(hruInfo$sub[isub] == rawOrdernum)
    smval = rawSM[, id + 1]
    smval[smval ==0] = NA
    smval[smval > 0.7] = NA # only take reliable sm range
    smvalAdj = smval * 50 # convert to mm unit
    smvalAdj[is.na(smvalAdj)] = -99
    sm4da[,hruInfo$lb[isub]:hruInfo$ub[isub]] = smvalAdj
    
  }
  
  
  datwu = matrix(rep(-99,nhru*nwu), ncol = nhru, nrow = nwu)
  sm4daAdj = rbind.data.frame(datwu, sm4da)
  

# read soil moisture from open loop simulation and plot temporal of SMAP and model--------

  
  olsm = read_fwf(file.path(maindir, folderName, '03openloop', 'Project','sm_1.dat'), skip = 0)
  olsm = data.frame(olsm)
  # location of 25th, 50th, and 75th of hru
  hruOrder = seq(1, nhru,1)
  lochru = c(round(quantile(hruOrder,0.25),0),
             round(quantile(hruOrder,0.50),0),
             round(quantile(hruOrder,0.75),0))
  ### PLOT before correction
  opFile25 = paste(figPath,'/',basinNames[ib],'_preC_25th_hru.pdf', sep = '')
  opFile50 = paste(figPath,'/',basinNames[ib],'_preC_50th_hru.pdf', sep = '')
  opFile75 = paste(figPath,'/',basinNames[ib],'_preC_75th_hru.pdf', sep = '')

  
  dat25 = data.frame(date = dateAll,
                     ol = olsm[,lochru[1]],
                     obs = sm4daAdj[,lochru[1]])
  dat50 = data.frame(date = dateAll,
                     ol = olsm[,lochru[2]],
                     obs = sm4daAdj[,lochru[2]])
  dat75 = data.frame(date = dateAll,
                     ol = olsm[,lochru[3]],
                     obs = sm4daAdj[,lochru[3]])
  # only use the simulation period
  dat25 = dat25[(nwu+1):nrow(dat25),]
  dat50 = dat50[(nwu+1):nrow(dat50),]
  dat75 = dat75[(nwu+1):nrow(dat75),]
  
  rownames(dat25) = NULL
  rownames(dat50) = NULL
  rownames(dat75) = NULL
  
  colObs = t_col("blue",percent = 55)
  colSim = t_col("red",percent = 55)
  
 
  # plot for 25th percentile
  dateVal = dat25$date
  uniqYear = unique(year(dat25$date))
  ny = length(uniqYear)
  
  pdf(opFile25,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat25$obs[idPlot]
    simPlot = dat25$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3,col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l",col=colSim) 
 
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
  # plot for 50th percentile
  dateVal = dat50$date
  uniqYear = unique(year(dat50$date))
  ny = length(uniqYear)
  
  pdf(opFile50,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat50$obs[idPlot]
    simPlot = dat50$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3,col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l",col=colSim) 
    
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
  # plot for 75th percentile
  dateVal = dat75$date
  uniqYear = unique(year(dat75$date))
  ny = length(uniqYear)
  
  pdf(opFile75,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat75$obs[idPlot]
    simPlot = dat75$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3,col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l",col=colSim) 
    
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
# correct smap data based on open loop -------------------------------------
  datol = data.frame(mm = month(dateAll),
                     dd = day(dateAll),
                     yy = year(dateAll),
                     olsm)
  sm4daAdj[sm4daAdj == -99] = NA
  datobs = data.frame(mm = month(dateAll),
                      dd = day(dateAll),
                      yy = year(dateAll),
                      sm4daAdj)
  
  summary(datol[,4], na.rm = T)
  summary(datobs[,4], na.rm = T)
  
  datolmean = daily_2_monthly_calc_mean(datol)
  datobsmean = daily_2_monthly_calc_mean(datobs)
  
  summary(datolmean[,1])
  summary(datobsmean[,1])
  
  datolsd = daily_2_monthly_calc_sd(datol)
  datobssd = daily_2_monthly_calc_sd(datobs)
  
  summary(datolsd[,1])
  summary(datobssd[,1])
  
  dateSpr = datobs[,1:3]
  monts = rep(seq(1,12,1),5) # create monthly data from 2015-2019
  
  datolssmean = monthly_seasonal_mean(datolmean, monts) 
  datobsssmean = monthly_seasonal_mean(datobsmean, monts)
  
  datolsssd = monthly_seasonal_mean(datolsd, monts)
  datobssssd = monthly_seasonal_mean(datobssd, monts)
  
  datolmeanAll = assing_value_monthly_2_daily(dateSpr, datolmean)
  datobsmeanAll = assing_value_monthly_2_daily(dateSpr, datobsmean)
  datolsdAll = assing_value_monthly_2_daily(dateSpr, datolsd)
  datobssdAll = assing_value_monthly_2_daily(dateSpr, datobssd)
  
  t = data.frame(datobs[,1:3], datobs[,4], datol[,4], 
                 datobsmeanAll[,1], datobssdAll[,1],
                 datolmeanAll[,1], datolsdAll[,1])
  colnames(t) = c('mm','dd','yy','obs','ol','obsmean','obssd','olmean','olsd')
  # normalization data
  sm4daAdj2 = datolmeanAll + (datobs[,-c(1,2,3)] - datobsmeanAll[,-c(1,2,3)])*(datolsdAll/datobssdAll)
  
  timeseriesAd = function(ts, tsol){
    ts[ts < 0 ] = 0
    maxval = max(tsol, na.rm = T)
    # assign Inf as NA
    ts[ts == Inf] = NA
    ts[ts > maxval] = maxval
    return(ts)
  }
  
  for(jj in 1:nhru){
    ts = sm4daAdj2[,jj]
    tsol = olsm[,jj]
    ts = timeseriesAd(ts, tsol)
    #summary(ts)
    sm4daAdj2[,jj] = ts
  }

  summary(olsm[,1])
  summary(sm4daAdj2[,1])
  ### PLOT after correction
  opFile25 = paste(figPath,'/',basinNames[ib],'_postC_25th_hru.pdf', sep = '')
  opFile50 = paste(figPath,'/',basinNames[ib],'_postC_50th_hru.pdf', sep = '')
  opFile75 = paste(figPath,'/',basinNames[ib],'_postC_75th_hru.pdf', sep = '')
  
  
  dat25 = data.frame(date = dateAll,
                     ol = olsm[,lochru[1]],
                     obs = sm4daAdj2[,lochru[1]])
  dat50 = data.frame(date = dateAll,
                     ol = olsm[,lochru[2]],
                     obs = sm4daAdj2[,lochru[2]])
  dat75 = data.frame(date = dateAll,
                     ol = olsm[,lochru[3]],
                     obs = sm4daAdj2[,lochru[3]])
  # only use the simulation period
  dat25 = dat25[(nwu+1):nrow(dat25),]
  dat50 = dat50[(nwu+1):nrow(dat50),]
  dat75 = dat75[(nwu+1):nrow(dat75),]
  
  rownames(dat25) = NULL
  rownames(dat50) = NULL
  rownames(dat75) = NULL
  
  colObs = t_col("blue",percent = 55)
  colSim = t_col("red",percent = 55)
  
  
  # plot for 25th percentile
  dateVal = dat25$date
  uniqYear = unique(year(dat25$date))
  ny = length(uniqYear)
  
  pdf(opFile25,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat25$obs[idPlot]
    simPlot = dat25$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3,col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l",col=colSim) 
    
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
  # plot for 50th percentile
  dateVal = dat50$date
  uniqYear = unique(year(dat50$date))
  ny = length(uniqYear)
  
  pdf(opFile50,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat50$obs[idPlot]
    simPlot = dat50$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3,col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l",col=colSim) 
    
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
  # plot for 75th percentile
  dateVal = dat75$date
  uniqYear = unique(year(dat75$date))
  ny = length(uniqYear)
  
  pdf(opFile75,width = 6,height = 6)
  par(mfrow=c(ny,1),mar=c(1.5,2,0.5,2),oma=c(0.1,1.5,0.1,1.5),mgp=c(2,0.5,0),tck=-0.03,cex.lab=0.9,cex.axis=0.8)
  for(p in 1:ny)
  {
    idPlot = which(year(dateVal)==uniqYear[p])
    obsPlot = dat75$obs[idPlot]
    simPlot = dat75$ol[idPlot]
    
    obsPlot[obsPlot == -99] = NA
    
    yrange= range(c(obsPlot,simPlot), na.rm = T)
    
    tickPos = seq.Date(as.Date(paste0(uniqYear[p],"-01-01")),as.Date(paste0(uniqYear[p],"-12-01")),"month")
    tickLable = format(tickPos,format = "%Y-%b")
    
    plot(dateVal[idPlot],obsPlot,type="p",pch = 3, col=colObs,ylim=yrange,
         xlab="",ylab="",xaxt="n",xaxs="i")
    points(dateVal[idPlot],simPlot,type="l", col=colSim) 
    
    
    axis(1,at=seq.Date(dateVal[1],dateVal[length(dateVal)],"months"),tck=-0.015,labels = F)
    axis(1,at=tickPos,tck=-0.03,labels = tickLable)
    axis(4)
    
    r2Dat = cor(simPlot,obsPlot,use = "na.or.complete")^2
    pbiasDat = pbias(simPlot,obsPlot)
    
    legend("topleft",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c(paste0("R-square: ",round(r2Dat,digits=2)),
                      paste0("pBias: ",round(pbiasDat,digits=2))),
           text.col = "darkred")
    
    legend("topright",bty="n",cex=0.7,text.font = 2,ncol=2,
           legend = c("Obs. SMAP","Openloop Sim"),lty=1,
           col = c(colObs,colSim))
  }
  mtext("soil moisture (mm)",2,outer=T,adj=0.5,cex=0.8)
  dev.off()
  
  
  # export bias correction SMAP data
  sm4daAdj2[is.na(sm4daAdj2)] = -99 # assign -99 for NA values
  par(mfrow = c(1,1))
  plot(seq(1, 1096,1),sm4daAdj2[nwu:1826,2], ylim = c(0,10), type = 'l')
  
  dqxhru = matrix(rep(0.08,(nsim+nwu)*nhru), ncol = nhru, nrow = (nsim+nwu)) # smap standard
  rfihru = matrix(rep(0,(nsim+nwu)*nhru), ncol = nhru, nrow = (nsim+nwu))
  
  #write.table(sm4daAdj2, file.path(maindir, folderName, '04dasta', 'Project','sm_obs_mn_std_monthly.txt'), 
  #            row.names = F, col.names = F)
  write.table(dqxhru, file.path(maindir, folderName, '04daens9','Project', 'dqx_hru.txt'), 
              row.names = F, col.names = F)
  write.table(rfihru, file.path(maindir, folderName, '04daens9', 'Project','rfi_hru.txt'), 
              row.names = F, col.names = F)
  
  write.table(sm4daAdj2, file.path(maindir, folderName, '04daens9', 'Project','sm_obs_mn_std_monthly.txt'), 
              row.names = F, col.names = F)
  
  #write.table(dqxhru, file.path(maindir, folderName, '05daens','Project', 'dqx_hru.txt'), 
  #            row.names = F, col.names = F)
  #write.table(rfihru, file.path(maindir, folderName, '05daens', 'Project','rfi_hru.txt'), 
   #           row.names = F, col.names = F)
  
}

