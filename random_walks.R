################################
## -------- R script ---------##
## ------ Random Walks ------ ##
################################

##step 1: preparations
  ##libraries needed
    library(gdistance)
    library(rgdal) 
  ##import rasters
DEM_300<-raster("D:\\Master_GIS\\R\\testraster\\iatrippa_tema_aggregated_fact10_arabmercator.tif")
LU_300<-raster("D:\\Master_GIS\\R\\testraster\\iatrippa_tema_LU_mod.tif")
  ##set start and goal locations
    start_name<-"Iatrippa"
    goal_name<-"Tema"
    A<-c(4401258,2791850)
    B<-c(4283116,3182311)
  ##cost function (inversing the algorithm to calculate conductance instead of resistance)
    herzog2012_wi <- function(s){1/(1337.8*s^6 + 278.19*s^5 - 517.39*s^4 - 78.199*s^3 + 93.419*s^2 + 19.825*s + 1.64)}
  ##number of neighbors
    N<-"4" 
  ##auxilliary function to calculate height difference between cells
    altdiff <- function(x){abs(x[2] - x[1])} 

## Step 2: calculate a transition matrix from the DEM

  ## 2a: crate an asymmetric transition matrix with the function applied to the raster, 
  ## for each cell to its N neighbours
  hd<-transition(DEM_300, altdiff, N, symm=T)
  
  ## 2b: calculate slope from heightdifferences
  slope<-geoCorrection(hd)
  
  ## 2c: create index for adjacent cells
  adj<-adjacent(DEM_300, cells=1:ncell(DEM_300), pairs=T, directions=N)
  
  ## 2d: calculate the energy expenditure of movement between adjacent cells using Herzog´s cost algorithm.
  speed<-slope
  speed[adj]<-herzog2012_wi(slope[adj])
  
  ## 2e: calculating conductance by relating speed with distance between cell centres
  ##     conductance is reciprocal of travel time
  conductance<-geoCorrection(speed)
  
  ##remove in-between-steps to free RAM
  rm(hd)
  rm(slope)
  rm(speed)
  rm(adj)
  gc()

## Step 3: correlate with LandUse raster

  ##create transition matrix from LandUse-raster using mean function
  ##take reciprocal to turn resistance into conductance
  tr_LU<-transition(1/LU_300, function(x) 1/mean(1/x), N, symm=T) ##using the harmonic mean
  tr_LU<-geoCorrection(tr_LU)
  tr_acc<-conductance*tr_LU

## step 4: calculate passage
  rw_A_B<-passage(tr_acc, A, B)
  rw_A_B_stan<-1-(log(1/rw_A_B)/cellStats(log(1/rw_A_B),max))
  ##invertieren um positive werte zu bekommen und logarithmieren um die gewichtung der kleinen werte zu erhöhen
  ##auf 0:1 skaliert indem man durch den max-wert dividiert
  ##durch das invertieren sind die werte nun umgedreht (höchste wahrsch bei 0, geringste bei 1), deswegen 1- am anfang um gegenwert zu erzeugen
  ##rw_A_B_stan2<-(asin(sqrt(rw_A_B_stan)))/cellStats(asin(sqrt(rw_A_B_stan)),max) ##arcsin transformation
  
  
  ##remove in-between-steps to free RAM
  rm(tr_LU)
  ##rm(tr_acc)
  gc()

  dir.create("results")
  writeRaster(rw_A_B, "./results/randomised_LCPs", format = "GTiff", overwrite=TRUE)
  writeRaster(rw_A_B_stan, "./results/randomised_LCPs_standartized", format = "GTiff", overwrite=TRUE)
