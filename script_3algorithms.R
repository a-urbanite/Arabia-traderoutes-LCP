##################################
##           R script           ##
##################################

##load necessary libraries
library(gdistance)
library(rgdal) 

## step 1:Import raster files and locations
  ##Raster
    DEM_300<-raster("D:\\Master_GIS\\R\\testraster\\iatrippa_tema_aggregated_fact10_arabmercator.tif")
    LU_300<-raster("D:\\Master_GIS\\R\\testraster\\iatrippa_tema_LU_mod.tif")
  ##Places
    A<-c(4401258,2791850) ##Iatrippa
    B<-c(4283116,3182311) ##Tema
  ##Cost functions
    ## Tobler1993 velocity
    tobler1993a  <-function(s){6*exp(-3.5*abs(s + 0.05))} # km/h
    tobler1993b <- function(s){0.36*exp(-3.5*abs(s + 0.05))} # m/min
    ## Minetti2002 metabolic costs J/(kg* m) walking
    minetti2002w    <- function(s){(280.5 *s^5 -58.7*s^4 - 76.8*s^3 + 51.9*s^2 + 19.6*s + 2.5)}
    minetti2002wi <- function(s){1/(280.5 *s^5 -58.7*s^4 - 76.8*s^3 + 51.9*s^2 + 19.6*s + 2.5)}
    ## Herzog2012 metabolic costs J/(kg*m) walking
    herzog2012_w    <- function(s){(1337.8*s^6 + 278.19*s^5 - 517.39*s^4 - 78.199*s^3 + 93.419*s^2 + 19.825*s + 1.64)}
    herzog2012_wi <- function(s){1/(1337.8*s^6 + 278.19*s^5 - 517.39*s^4 - 78.199*s^3 + 93.419*s^2 + 19.825*s + 1.64)}
    ##Bell and Lock
    ##BellLock_i <- function(s){tan(s)/tan(1°)}
  ##auxilliary funciton to calculate heigth difference between cells
    altdiff <- function(x){abs(x[2] - x[1])} 
  ##Number of Neighbors
    N="8"

## Step 2: calculations independent from cost function
  ##calculate slope
    slope<-geoCorrection(transition(DEM_300, altdiff, N, symm=F))
  ##set index for adjacent cells
    adj<-adjacent(DEM_300, cells=1:ncell(DEM_300), pairs=T, directions=N)
  ##calculate the transition matrix of the LandUse-raster 
    tr_LU<-transition(1/LU_300, function(x) 1/mean(1/x), N, symm=F) ## taking the harmonic mean instead of the arithmetic one as proposed by van Etten
    tr_LU<-geoCorrection(tr_LU) 
    
## Step 3: calculate a transition matrix from the DEM, correlate it with the Landuse Raster adn calculate the routes
  ## 2a Tobler Hiking function
    speed<-slope
    speed[adj]<-tobler1993a(slope[adj])    ##apply Tobler function to [adj] selection of cells of slope
    tr_DEM<-geoCorrection(speed)           
    tr_acc<-tr_DEM*tr_LU                   ##correlate with Landuse-Raster
    rm(tr_DEM)
    gc()
    ##calculate the LCP´s
    Tobler_A_to_B<-shortestPath(tr_acc, A, B, output="SpatialLines")
    Tobler_B_to_A<-shortestPath(tr_acc, B, A, output="SpatialLines")
    ##free memory
    rm(tr_acc)
    gc()

  ## 2b Minetti2002
    speed<-slope
    speed[adj]<-minetti2002wi(slope[adj])
    tr_DEM<-geoCorrection(speed)
    tr_acc<-tr_DEM*tr_LU
    
    rm(tr_DEM)
    gc()
  
    Minetti_A_to_B<-shortestPath(tr_acc, A, B, output="SpatialLines")
    Minetti_B_to_A<-shortestPath(tr_acc, B, A, output="SpatialLines")
  
    rm(tr_acc)
    gc()

  ## 2c Herzog2012
    speed<-slope
    speed[adj]<-herzog2012_wi(slope[adj])
    tr_DEM<-geoCorrection(speed)
    tr_acc<-tr_DEM*tr_LU
    
    rm(tr_DEM)
    gc()
    
    Herzog_A_to_B<-shortestPath(tr_acc, A, B, output="SpatialLines")
    Herzog_B_to_A<-shortestPath(tr_acc, B, A, output="SpatialLines")
    
    rm(tr_LU)
    gc()

## step 4: plotting, clear memory and export

  plot(raster(tr_acc), xlab="x coordinate(m)", ylab="y coordinate(m)", legend.lab="conductance", main="method 3c")
  lines(Herzog_A_to_B, col="red", lwd=2)
  lines(Herzog_B_to_A, col="blue")
  lines(Minetti_A_to_B, col="red", lwd=2)
  lines(Minetti_B_to_A, col="blue")
  lines(Tobler_A_to_B, col="red", lwd=2)
  lines(Tobler_B_to_A, col="blue")
  
  rm(tr_acc)
  gc()
  
  dir.create("results")
  writeOGR(SpatialLinesDataFrame(Tobler_A_to_B, data.frame(id=1:length(Tobler_A_to_B))), "./results", "Tobler_A_to_B", driver="ESRI Shapefile")
  writeOGR(SpatialLinesDataFrame(Tobler_B_to_A, data.frame(id=1:length(Tobler_B_to_A))), "./results", "Tobler_B_to_A", driver="ESRI Shapefile")
  writeOGR(SpatialLinesDataFrame(Minetti_A_to_B, data.frame(id=1:length(Minetti_A_to_B))), "./results", "Minetti_A_to_B", driver="ESRI Shapefile")
  writeOGR(SpatialLinesDataFrame(Minetti_B_to_A, data.frame(id=1:length(Minetti_B_to_A))), "./results", "Minetti_B_to_A", driver="ESRI Shapefile")
  writeOGR(SpatialLinesDataFrame(Herzog_A_to_B, data.frame(id=1:length(Herzog_A_to_B))), "./results", "Herzog_A_to_B", driver="ESRI Shapefile")
  writeOGR(SpatialLinesDataFrame(Herzog_B_to_A, data.frame(id=1:length(Herzog_B_to_A))), "./results", "herzog_B_to_A", driver="ESRI Shapefile")
