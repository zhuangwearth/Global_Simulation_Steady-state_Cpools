####################################################################
###    Global Soil Carbon Simulation with MIMICSa and MIMICSc   ####
####################################################################

library(raster)
library(oceanmap)
library(rootSolve)
library(FME)

#################################################################
###                       load NPP and soil data            #####
#################################################################
setwd('D:/Csoil/steay_simu/gC_per_kg_1m/')
# The model has five driving variables
# V1: ANPP
# V2: TSOI
# V3: CLAY  
# V4: bulkD  
# V5: Qmax

### Load V1 NPP -----
clmpath = 'D:/Data/Csoil/steady_simu/gC_per_kg_1m/inputData/dataCLM'
## NPP is 0.5 degree
npp_map = raster(paste0(clmpath, 'npp_gC_per_m2_y.nc4'))
# divide 2 to get ANPP
npp_map = npp_map/2 
# transform raster to matrix
npp_matrix = raster2matrix(npp_map)

### load V2 soil temperature data   -----
## Tsoil is 0.5 degree
tsoil_map = raster(paste0(clmpath, 'SoilTemp1m.nc'))
# transform raster to matrix
tsoil_matrix = raster2matrix(tsoil_map)

### load V3 CLAY content (%)   -----
hwsdpath = 'D:/Data/Csoil/steady_simu/gC_per_kg_1m/inputData/HWSD'
# Please convert T_CLAY_HWSD.tif from 0.05 to 0.5 before load
CLAY_map   <- raster(paste0(hwsdpath, 'CLAY_1m_0.5deg.tif'))
# NOTE! quality control is needed because CLAY HAS ZERO VALUES
# transform raster to matrix
CLAY_matrix = raster2matrix(CLAY_map)

### load V4 bulk density  -----
bulkD_map   <- raster(paste0(hwsdpath, 'bulkD_1m_0.5deg.tif'))
# transform raster to matrix
bulkD_matrix = raster2matrix(bulkD_map)

### load V5 Qmax -----
# MIMICSc is different from MIMICSb HERE.
# MOCmax based on observation is for Qmax for MIMICSc. 
pathmoc = 'D:\\Csoil\\MAOM_Globe\\rawdata\\mocMax1m_mg_per_cm3.nc'
Qmax_map   <- raster(pathmoc)
# transform raster to matrix
Qmax_matrix = raster2matrix(Qmax_map)

#################################################################
###       conduct Csoil simulation with matrix data         #####
#################################################################
# simulate at one grid
ilon = 91.75	# 45.75 
ilat = 37.75  # 20.25

ilon = 179.75 
ilat = 71.25 
ix = (ilon + 180.25)*2
iy = (ilat + 90.25)*2
ix;
iy
ANPP   <- npp_matrix[ix, iy]      # gC per m2 per year. typical range: 110 ~ 1500
TSOI   <- tsoil_matrix[ix, iy]    # soil temperature
CLAY   <- CLAY_matrix[ix, iy]     # CLAY
bulkD  <- bulkD_matrix[ix, iy]    # bulk density (g cm-3)
Qmax   <- Qmax_matrix[ix, iy]  

depth  <- 100                     # depth of soil carbon 

source('code/cal_steady_MIMICSab.R')
print(beta)

res_tmp <- matrix(nrow = 2,ncol = 7)

test <- tryCatch({
  steady_func(ANPP,TSOI,CLAY,bulkD, Qmax)
}, error = function(e) res_tmp)

test[[1]]
test[[2]]
sum(test[[1]])
sum(test[[2]])
remove(steady_func)


# define function to introduce the matrix to models
func_test <- function(ix,iy,npp_matrix,tsoil_matrix,CLAY_matrix,bulkD_matrix, Qmax_matrix){
  ANPP  <- npp_matrix[ix, iy]      # gC per m2 per year. typical range: 110 ~ 1500
  TSOI   <- tsoil_matrix[ix, iy]    # soil temperature
  CLAY   <- CLAY_matrix[ix, iy] 
  bulkD  <- bulkD_matrix[ix, iy]
  Qmax   <- Qmax_matrix[ix, iy] 
  
  res <- matrix(nrow = 2,ncol = 7)
  # In this case, C steady have no value
  if (ANPP = 0 | CLAY == 0 | any(is.na(c(ANPP,TSOI,CLAY,bulkD)))) 
    {sa = rep(NA, 7)
    sb = rep(NA, 7)
    }
  # In this case, C steady have values
  else
    {source('./code/cal_steady_MIMICSab.R')
      res_ab <- tryCatch({
                          steady_func(ANPP,TSOI,CLAY,bulkD, Qmax)
                         }, error = function(e) res_tmp)
      sa = res_ab[[1]]
      sb = res_ab[[2]]
      remove(steady_func)
    }
  res[1,] <- sa
  res[2,] <- sb
  return(res)
}

locs <- matrix(nrow = 720*360,ncol = 2)
i <- 1
for(ix in seq(1:720)){
  for(iy in seq(1:360)){
    locs[i,] <- c(ix,iy)
    i=i+1
  }
}

### conduct calculation
library(doParallel)
library(foreach)

detectCores()

myCluster <- makeCluster(60) # type of cluster
# calculate
registerDoParallel(myCluster)

system.time(
  output <-foreach(i=1:259200,.combine = rbind)%dopar%{
    ix <- locs[i,][1]
    iy <- locs[i,][2]
    ilon = ix/2 - 180.25
    ilat = iy/2 - 90.25
    print(paste0(ilon, ' + ', ilat))
    res <- func_test(ix,iy,npp_matrix,tsoil_matrix,CLAY_matrix,bulkD_matrix, Qmax_matrix)
    # sa <- res[1,]
    # sb <- res[2,]
   
    c(ilon,ilat,res[1,],res[2,])
  }
)

write.csv(output,'output/Cpools_MIMICSa_MIMICSc.csv', row.names = F)
stopCluster(myCluster)
### End calculation

# combine results of the three versions of MIMICS model
MaMb = read.csv('output/Cpools_MIMICSa_MIMICSb.csv')
MaMbMc = rbind(MaMb[,1:16], output[,10:16])
write.csv(MaMbMc,'output/Cpools_MIMICSabc.csv', row.names = F)