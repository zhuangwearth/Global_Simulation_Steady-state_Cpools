# The code file is called by the main files "simulate_global_Csoil_MIMICSa_MIMICSb.R"
# and 'simulate_global_Csoil_MIMICSc_MOCmax.R'

###
###  Solve steady state for a system of non-linear equations in MIMICSab model
###
# wei 14/Jun/2022
# revised according to ode function in Rose Abramoff 5/4/2017 and 
# orginal codes from Wieder et al. Geosci. Model Dev. Discuss., 8, 2011 Created Aug 16 2013 (Modified May 2014)
# Here tao modified by site level productivity (sqrt(ANPP/100))

# Single layer model with: 
# 2 litter C pools (LIT), corresponding to metabolic and structural litter.
# 2 microbial pool (MIC; i.e. R vs. K strategists)
# 3 SOM pools corresponding to physically & chemically protected & available pools
# ADDS COMPLEXITY TO LITTER DECOMPOSITION

### Three processing steps 
# (0) Prepare parameter and model
# (A) Reads in site level data  the specific sites
# (B) Calculates steady state C pools using XEQ & STODE function

####################################################################### 
###               (0) Prepare parameter and model                   ###
#######################################################################
source("./code/model_parameter_prepare.R")

####################################################################### 
###      (A) Reads in site level data from the specific sites       ###
####################################################################### 
### Driving variables of steady state function
# ANPP
# TSOI
# CLAY

steady_func <- function(ANPP,TSOI,CLAY,blukD, Qmax){
  
  source("./code/model_parameter_prepare.R")
  # variables from grid of global data
  # ANPP   <- ANPP            # gC per m2 per year. typical range: 110 ~ 1500
  # TSOI   <- TSOI            # soil temperature, adopted from drivering data 
  # CLAY   <- CLAY # CLAY
  # bulkD  =  bulkD              # bulk density (g cm-3)
  
  
  ## EST_LIT is related with NPP
  EST_LIT_in  <- ANPP / (365*24)   		# input litter amount, to gC/m2/h  from gC/m2/y
  
  Tao_MOD1 <- sqrt(ANPP/100)              # basicaily standardize against NWT 
  
  # limit the range of Tao_MOD1 in 0.8 ~ 1.2
  if (Tao_MOD1 < 0.8)
     Tao_MOD1 = 0.8

  if (Tao_MOD1 > 1.2)
     Tao_MOD1 = 1.2  
  
  
  # Qmax in MIMICSb
  # Qmax  = 10^(0.297 * log10(CLAY) + 3.355)/1000 # unit, g C kg-1 Soil)
  # Qmax  = bulkD * Qmax                          # mg C per cm3 Soil
  
  ####################################################################### 
  ### (B) Calculates steady state C pools using XEQ & STODE function  ###
  ####################################################################### 
  #
  # initialize pools
  int.val = c(1, 1)
  LIT     <- int.val 
  MIC     <- int.val/100
  SOM     <- rep(NA, 3) 
  SOM[1]  <- int.val[1]
  SOM[2]  <- int.val[2]
  SOM[3]  <- int.val[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  # initial states             
  Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
              MIC_1 = MIC[1], MIC_2 = MIC[2], 
              SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
  
  ####################################################################### 
  ### (B.1) Calculates steady state for MIMICSa                       ###
  ####################################################################### 
  #
  ## calculate the values of mid-parameters
  param_mid_MIMICSa = cal_pars_MIMICSa(EST_LIT_in, depth, TSOI, CLAY, 
                                       LIG, CN, 
                                       Vslope, Vint, aV, MOD1, Kslope, Kint, aK, Tao_MOD1)
  attach(param_mid_MIMICSa)
  
  # parameters
  Tpars_MIMICSa   <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                        fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                        tao = tao, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                        desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  
  # solve steady state for MIMICSa
  Csteady_MIMICSa  <- stode(y = Ty, time = 1e6, fun = XEQ_MIMICSa, parms = Tpars_MIMICSa, positive = TRUE)
  
  
  ####################################################################### 
  ### (B.2) Calculates steady state for MIMICSb                       ###
  ####################################################################### 
  #
  ## calculate the values of mid-parameters
  param_mid_MIMICSb = cal_pars_MIMICSb(EST_LIT_in, depth, TSOI, CLAY, 
                                       LIG, CN, 
                                       Vslope, Vint, aV, MOD1, Kslope, Kint, aK, Tao_MOD1, Qmax)
  attach(param_mid_MIMICSb)
  
  # parameters
  Tpars_MIMICSb   <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
                        fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
                        tao = tao, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
                        desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO, Qmax=Qmax)
  
  # solve steady state for MIMICSb
  Csteady_MIMICSb  <- stode(y = Ty, time = 1e6, fun = XEQ_MIMICSb, parms = Tpars_MIMICSb, positive = TRUE)
  
  # Store the final varible that we need:
  # the C pools at steady state
  # Csteady_MIMICSa[[1]]
  # Csteady_MIMICSb[[1]]
  return(c(Csteady_MIMICSa,Csteady_MIMICSb))
}