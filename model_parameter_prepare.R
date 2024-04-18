# The code file is called by the main files "simulate_global_Csoil_MIMICSa_MIMICSb.R"
# and 'simulate_global_Csoil_MIMICSc_MOCmax.R'

library(rootSolve)
library(FME)


####################################################################### 
###               (0) Prepare parameter and model                   ###
#######################################################################

## (0.1) parameters or constants in mimicsa
depth  <- 100 # depth_cm    # depth of soil carbon

beta     <- 1.05
FI       <- c(0.05, 0.05)	  
Vslope   <- array(0.063,dim=6)
Vint     <- 5.47
aV       <- 8e-6
MOD1     <- c(10, 2, 10, 3, 3, 2) 

Kslope   <- array(NA,dim=6)
Kslope[1]<- 0.017 #META LIT to MIC_1
Kslope[2]<- 0.027 #STRU LIT to MIC_1 
Kslope[3]<- 0.017 #AVAI SOM to MIC_1 
Kslope[4]<- 0.017 #META LIT to MIC_2
Kslope[5]<- 0.027 #STRU LIT to MIC_2
Kslope[6]<- 0.017 #AVAI SOM to MIC_2
Kint     <- 3.19
aK       <- 10
KO       <- c(4,4)      #scalar modifies Km of Oxidat	  
# CUE
CUE      <- c(0.55, 0.25, 0.75, 0.35)  #for LITm and LITs entering MICr and MICK, respectively

times_fc = 1

## parameter CN
LIG    <- 21              
CN     <- 49

######### changed k and a, to wieder 2015 A3
k        <- 2    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
a        <- 2    #2			#increased from 4.0 to 4.5

## (0.2)  sub-model to calculate the 25 mid-parameters in mimics
# construct a sub-model to calculate the 25 mid-parameter directly used in the pool change functions
# model input: basic parameter
# model output: mid-parameter for XEQ
# note: The units of I, VMAX, tao and desorb are converted from per hour to per year in the function
  
cal_pars_MIMICSa <- function(EST_LIT_in, depth, TSOI, CLAY, 
                     LIG, CN, 
                     Vslope, Vint, aV, MOD1, Kslope, Kint, aK, Tao_MOD1) {	
   	fCLAY  <- CLAY/100			 # clay content, convert from clay fraction to %
						 		
    # Calculate I[1:2], FI[1:2] parameters related to litter chane
    # I=f(x=c(ANPP, depth, LIG, CN, CLAY), pars= c( ) )
    calCN  <- (1 / CN) / 2.5 * 100 
    fMET   <- 0.85 - 0.013 * LIG/calCN 			#as calculated in DAYCENT
    
    EST_LIT  <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
    Inp        <- rep(NA, 2)
    Inp[1]     <- (EST_LIT / depth) * fMET   #partitioned to layers
    Inp[2]     <- (EST_LIT / depth) * (1-fMET)
    Inp        <- Inp * 24 * 365              # convert unit from per hour to per year
    
    # Calculate Vmax(T) & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
    # VMAX = f(x = TSOI, pars = c(Vslope, Vint, aV, MOD1))
    Vmax     <- exp(TSOI * Vslope + Vint) * aV
    VMAX     <- Vmax * MOD1  # !!MODIFIERS AS IN MIMICS2_b!!
    VMAX     <- VMAX * 24 * 365    # convert unit from per hour to per year
    
    #Calculate Km(T) 	
    # Km = f(x = TSOI, pars = c(Kslope, Kint, aK, MOD2, pSCALAR)) 
    Km       <- exp(Kslope * TSOI + Kint) * aK
	
	######### changed k and a, to wieder 2015 A3
    # k        <- 3    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
    # a        <- 0.8    #2.2			#increased from 4.0 to 4.5
    pSCALAR  <- 1.0/ (a * exp(-k*(sqrt(fCLAY)))) #Scalar for texture effects on SOMp
    MOD2     <- c( 0.125, 0.5 , 0.25 * pSCALAR, 0.5, 0.25, pSCALAR/6)   
    KM       <- Km * MOD2	
    
    # parameter related to mic turnover
    # tao = f(x = c(LIG, CN, fCLAY, ANPP)   
    tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	# microbial biomass turover rates of r and K types
    tao      <- tao * Tao_MOD1    
    tao      <- tao * 24 * 365    # convert unit from per hour to per year
    
    fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
    fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
	
	######### changed 4 times to wieder 2015 A3
	fCHEM    <- times_fc * fCHEM  
    fAVAI    <- 1- (fPHYS + fCHEM)
    
    # desorb
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!! 
    desorb   <- desorb * 24 * 365   # convert unit from per hour to per year
      
    return(list(Inp=Inp, VMAX=VMAX, KM=KM, KO=KO, tao=tao, 
                fPHYS=fPHYS, fCHEM=fCHEM, fAVAI=fAVAI, desorb=desorb)) 
    }
# End of sub-model to calculate the 25 mid-parameters


cal_pars_MIMICSb <- function(EST_LIT_in, depth, TSOI, CLAY, 
                     LIG, CN, 
                     Vslope, Vint, aV, MOD1, Kslope, Kint, aK, Tao_MOD1, Qmax) {	
   	fCLAY  <- CLAY/100			 # clay content, convert from clay fraction to %
						 		
    # Calculate I[1:2], FI[1:2] parameters related to litter chane
    # I=f(x=c(ANPP, depth, LIG, CN, CLAY), pars= c( ) )
    calCN  <- (1 / CN) / 2.5 * 100 
    fMET   <- 0.85 - 0.013 * LIG/calCN 			#as calculated in DAYCENT
    
    EST_LIT  <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
    Inp        <- rep(NA, 2)
    Inp[1]     <- (EST_LIT / depth) * fMET   #partitioned to layers
    Inp[2]     <- (EST_LIT / depth) * (1-fMET)
    Inp        <- Inp * 24 * 365              # convert unit from per hour to per year
    
    # Calculate Vmax(T) & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
    # VMAX = f(x = TSOI, pars = c(Vslope, Vint, aV, MOD1))
    Vmax     <- exp(TSOI * Vslope + Vint) * aV
    VMAX     <- Vmax * MOD1  # !!MODIFIERS AS IN MIMICS2_b!!
    VMAX     <- VMAX * 24 * 365    # convert unit from per hour to per year
    
    #Calculate Km(T) 	
    # Km = f(x = TSOI, pars = c(Kslope, Kint, aK, MOD2, pSCALAR)) 
    Km       <- exp(Kslope * TSOI + Kint) * aK
	
	######### changed k and a, to wieder 2015 A3
    # k        <- 3      #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
    # a        <- 0.8    #2.0			#increased from 4.0 to 4.5
    pSCALAR  <- 1.0/ (a * exp(-k*(sqrt(fCLAY)))) #Scalar for texture effects on SOMp
    MOD2     <- c( 0.125, 0.5 , 0.25 * pSCALAR, 0.5, 0.25, pSCALAR/6)   
    KM       <- Km * MOD2	
    
    # parameter related to mic turnover
    # tao = f(x = c(LIG, CN, fCLAY, ANPP)   
    tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	# microbial biomass turover rates of r and K types
    tao      <- tao * Tao_MOD1    
    tao      <- tao * 24 * 365    # convert unit from per hour to per year
    
    fPHYS    <- c(0, 0)   	                                    #fraction to SOMp
    fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
	
	######### 1 times 
	fCHEM    <- times_fc * fCHEM  
    fAVAI    <- 1- (fPHYS + fCHEM)
    
    # desorb
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))      #CHANGED FOR GLOBAL RUN!!! 
    desorb   <- desorb * 24 * 365   # convert unit from per hour to per year
    
    Qmax     <- Qmax
	
    return(list(Inp=Inp, VMAX=VMAX, KM=KM, KO=KO, tao=tao, 
                fPHYS=fPHYS, fCHEM=fCHEM, fAVAI=fAVAI, desorb=desorb, Qmax = Qmax)) 
    }
# End of sub-model to calculate the 25 mid-parameters


## (0.3) Function used to calculate steady state C pools using XEQ & STODE function

# Set up model function for MIMICSa
XEQ_MIMICSa <- function(t, y, pars) { #t = time, y = state, pars = pars
	with (as.list(c(y, pars)),{
		# define parameters
		# basic parameters need sensitivity analysis
					
		#Flows to and from MIC_1
		LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
		LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
		
		tMIC_1    = MIC_1^beta * tao[1] 
		MICtrn[1] = tMIC_1 * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
		MICtrn[2] = tMIC_1 * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
		MICtrn[3] = tMIC_1 * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
 
		SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1

		#Flows to and from MIC_2
		LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #MIC_2 decomp of MET litter
		LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #MIC_2 decomp of SRUCTURAL litter
		
		tMIC_2    = MIC_2^beta * tao[2] 
		MICtrn[4] = tMIC_2 * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
		MICtrn[5] = tMIC_2 * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
		MICtrn[6] = tMIC_2 * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  

		SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
		
		DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
		OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
					  (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
			#can make fluxes from CHEM a function of microbial biomass size?

		dLIT_1 = Inp[1] - Inp[1]*FI[1] - LITmin[1] - LITmin[3]                              # Met Litter output: SOM_1, mic_r_K use
		dLIT_2 = Inp[2] - Inp[2]*FI[2] - LITmin[2] - LITmin[4]                              # Structral Litter output: SOM_2, mic_r_K use

		dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])  # mic_r
		dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  # mic_k 

		dSOM_1 = Inp[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb                             # PHYS protected SOM 			
		dSOM_2 = Inp[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT                            # chemically lazy som
		dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]       # available SOM

		list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
		})
}


# Set up model function for MIMICSb
XEQ_MIMICSb <- function(t, y, pars) { #t = time, y = state, pars = pars
	with (as.list(c(y, pars)),{
				
		#Flows to and from MIC_1
		LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
		LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
        tMIC_1    = MIC_1^beta * tao[1] 
		MICtrn[1] = tMIC_1 * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
		MICtrn[2] = tMIC_1 * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
		MICtrn[3] = tMIC_1 * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
		SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1

		#Flows to and from MIC_2
		LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #MIC_2 decomp of MET litter
		LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #MIC_2 decomp of SRUCTURAL litter
        tMIC_2    = MIC_2^beta * tao[2] 
		MICtrn[4] = tMIC_2 * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
		MICtrn[5] = tMIC_2 * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
		MICtrn[6] = tMIC_2 * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
		SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
		
		### move to 'simulate_global_Csoil...R'
		### revised
		# Qmax in MIMICSb
		# Qmax  = 10^(0.297 * log10(CLAY) + 3.355)/1000 # unit, g C kg-1 Soil)
		# Qmax  = bulkD * Qmax  
		
        DEsorb = SOM_1 * desorb                		#desorbtion of PHYS to AVAIL (function of fCLAY)
        
		ADsorb = SOM_3 * 6 * desorb * (1 - SOM_1 / Qmax)    # different from MIMICSa
        
        OXIDAT = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
					  (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
			#can make fluxes from CHEM a function of microbial biomass size?

		dLIT_1 = Inp[1] - Inp[1]*FI[1] - LITmin[1] - LITmin[3]                              # Met Litter output: SOM_1, mic_r_K use
		dLIT_2 = Inp[2] - Inp[2]*FI[2] - LITmin[2] - LITmin[4]                              # Structral Litter output: SOM_2, mic_r_K use

		dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])  # mic_r
		dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  # mic_k 
        
        # revised
		dSOM_1 = Inp[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb +  ADsorb                            # PHYS protected SOM 			
		dSOM_2 = Inp[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT                            # chemically lazy som
		dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2] - ADsorb      # available SOM

		list(c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3))
		})
}
