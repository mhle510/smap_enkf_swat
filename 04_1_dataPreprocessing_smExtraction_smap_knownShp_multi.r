####################################################################################
# Manh-Hung Le - 2021 June 29
# objectives:
####################################################################################
library(lubridate)
library(tidyverse)
library(raster)
library(sp)
library(rgdal)
library(readxl)
library(tidyverse)
library(lubridate)
library(rgeos)
library(tictoc)
library(parallel)
library(doParallel)
#library(sf)

# symbol - full name
# chu - Chu
# tso - ThanhSon
# vye - VinhYen
# kle - KheLech
# mcc - Mucangchai
# nhu - NaHu
# qcu - Quang Cu
# lso - Lang Son
# bye - Ban Yen
# gso - Giang Son
# bli - Binh Lieu
# tnh - Thuong Nhat
# gvo - Gia Vong
# ach - An Chi
# aho - An Hoa

# swat folder
maindir = 'D:/sda'
basinNames = c('gvo','aho','bye','slu', 'chu','gso','nkh','xla')
riverSystems = c('sk', 'sk','mk','sk','ht','mk','ca','ma')
nb = length(basinNames)

# smap pm folder
smtifPath = 'E:/GlobalSM/processedVnbasins/9km/pm'
smFiles = list.files(smtifPath, pattern = '.tif', full.names = T,recursive = T)
# select 2017 data onward
date = basename(smFiles) %>% str_remove('.tif')
locid = which(substr(date, 1, 4) == '2017')
locidlast = which(substr(date, 1, 4) == '2020')


smFilesAdj = smFiles[locid[1]:(locidlast[1]-1)]
dateAdj = date[locid[1]:(locidlast[1]-1)]

smFilesAdj[1]
dateAdj[1]

smFilesAdj[length(smFilesAdj)]
dateAdj[length(smFilesAdj)]

nd = length(smFilesAdj)

dateF = make_date(year = substr(dateAdj,1,4), month = substr(dateAdj,5,6), day = substr(dateAdj,7,8))

# own-built packages
codePath = 'D:/sda/0_inputPreparation/rCode'
source('C:/Users/Admin/Dropbox/RStudy/0.Code/built_in/sm_support_functions.r')

for(ib in c(1,3,5)){
  tic()
  mainPath = paste(maindir,'/','swat_',riverSystems[ib],'_',basinNames[ib], sep = '')
  catname = basinNames[ib]
  cat('========= working with catchment ==========','\n')
  cat(catname,'\n')
  basinShp = readOGR(paste(mainPath,'/','00dataPreparation','/','knownSubbasins','/',
                           catname,'_subs.shp', sep = ''))
  #basinSf = read_sf(paste(mainPath,'/','00dataPreparation','/','knownSubbasins','/',
  #                       catname,'_subs.shp', sep = ''))
  # create new folder to store precipitation at centroid sub-basins
  dir.create(file.path(mainPath, '00dataPreparation','processedSM'), showWarnings = F)
  processedPath = file.path(mainPath, '00dataPreparation','processedSM')
  
  CRS.wgs = CRS("+init=epsg:4326")
  basinShpwgs = spTransform(basinShp, CRS.wgs)
  nameSub = basinShpwgs$Subbasin
  nShp = length(basinShpwgs)
  gridPoints = data.frame(name = paste('sm',
                                       paste0(formatC(as.numeric(nameSub),width = 3,flag = 0)),sep = ''),
                          X = getSpPPolygonsLabptSlots(basinShpwgs)[,1],
                          Y = getSpPPolygonsLabptSlots(basinShpwgs)[,2])
  #cat('========= working with subbasin ==========','\n')
  smDat = mat.or.vec(nd, nShp)
  
  #noCores = detectCores() - 12
  #cl1 = makeCluster(noCores)
  
  no_cores = detectCores(logical = TRUE)
  cl = makeCluster(13)  
  #registerDoParallel(cl)  
  
  # export library raster to all cores
  clusterEvalQ(cl, library("raster"))
  clusterEvalQ(cl, library("rgeos"))
  # convert basin from SpatialPolygonsDataFrame to SpatialPolygons (required format by the extract() function)
  basinGeom = geometry(basinShpwgs)
  
  basinShpwgsSingle = basinShpwgs
  # assign 1 for polygonID
  basinShpwgsSingle$PolygonId = 1
  # dissovle automatically basinShpwgsSingle
  rgeos::set_RGEOS_CheckValidity(2L)
  basinShpwgsSingle = gUnaryUnion(basinShpwgsSingle, id = basinShpwgsSingle$PolygonId)
  plot(basinShpwgsSingle)
  tfrac = parLapply(cl, 1:nd,
                     function(ll, smFilesAdj, basinShpwgsSingle, calc_fraction_area){
                       rawRas = raster(smFilesAdj[[ll]]) # read raster in this step
                       #plot(rawRas, xlim = c(106.6, 106.99), ylim = c(16.7, 16.99))
                       #plot(basinShpwgsSingle, add = T)
                       ratioDat = calc_fraction_area(rawRas, rawShp = basinShpwgsSingle)
                      }, smFilesAdj, basinShpwgsSingle, calc_fraction_area
  )
  smFrac = do.call(rbind,tfrac)
  
  tic()
  t = parLapply(cl = cl, 1:nd,
                function(ll,smFilesAdj,basinGeom, calc_weightMatrix_sm){
                  # note that we need to explicitly 
                  #declare all "external" variables used inside clusters
                  rawRas = raster(smFilesAdj[[ll]]) # read raster in this step
                  opDat = lapply(1:length(basinGeom),FUN = function(i){
                    calc_weightMatrix_sm(rawShp = basinGeom[i,], rawRas)
                  })
                  # merge data from list to a vector
                  do.call(c,opDat)
                }, smFilesAdj,basinGeom, calc_weightMatrix_sm# you need to pass these variables into the clusters
  )
  toc()
  
  smDat = do.call(rbind,t) # each subcatchment is a column
  
  smDat = data.frame(date = dateF,
                      smDat)
  
  colnames(smDat) = c('date',
                       paste('sm',paste0(formatC(as.numeric(nameSub),width = 3,flag = 0)),sep = ''))
  pointsCoord = gridPoints[,c(2:3)]
  
  write.csv(smFrac, paste(processedPath,'/','smFrac.csv', sep = ''), row.names = F)
  write.csv(smDat, paste(processedPath,'/','smSMAP.csv', sep = ''), row.names = F)
  write.csv(pointsCoord, paste(processedPath,'/','smSMAPcoord.csv', sep = ''), row.names = F)
  toc()
}
