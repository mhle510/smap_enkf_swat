# RUN AFTER RUN OPEN LOOP !!!!!!!

#install.packages('reader')
library(reader)
library(rgdal)
library(tidyverse)

# swat folder
maindir = 'D:/sda'
basinNames = c('vye','gvo','aho','bye','slu', 'chu','gso','nkh','xla')
riverSystems = c('ht','sk', 'sk','mk','sk','ht','mk','ca','ma')
# create folder

nb = length(basinNames)

for(ib in c(1,3,4,5,8)){
  folderName = paste('swat','_',riverSystems[ib],'_',
                     basinNames[ib], sep = '')
  basinShp = readOGR(paste(maindir,'/',folderName,'/','00dataPreparation','/','knownSubbasins','/',
                           basinNames[ib],'_subs.shp', sep = ''))
  nSubs = length(basinShp)
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
    smvalAdj = smval * 5 # convert to mm unit
    smvalAdj[is.na(smvalAdj)] = -99
    sm4da[,hruInfo$lb[isub]:hruInfo$ub[isub]] = smvalAdj
    
  }
  
  
  dqxhru = matrix(rep(0.04,nsim*nhru), ncol = nhru, nrow = nsim) # smap standard
  rfihru = matrix(rep(0,nsim*nhru), ncol = nhru, nrow = nsim)
  
  write.table(sm4da, file.path(maindir, folderName, '03darun', 'sm_obs_mn_std_monthly.txt'), 
              row.names = F, col.names = F)
  write.table(dqxhru, file.path(maindir, folderName, '03darun', 'dqx_hru.txt'), 
              row.names = F, col.names = F)
  write.table(rfihru, file.path(maindir, folderName, '03darun', 'rfi_hru.txt'), 
              row.names = F, col.names = F)
}
