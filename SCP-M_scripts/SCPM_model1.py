#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
#The [simple carbon project] model, in process
# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

# Working directory for results
ResultsDir='Results/'
LastDir='Last_run/'
BatchResultsDir='Results/Batch/'


#Import modules
import math
import numpy as np
import matplotlib.pyplot as plt
import os


#Import matrices and parameters
from SCPM_matrices import (makeBiomat_out,makeBiomat_in,
makeZarr,makeT1mat, makeT2mat,makeE1mat, makeE2mat, makeFXarr)
from SCPM_parameters import (Psi1,Psi2,gamma1,gamma2,Zed,bScale,Fract,PV0,PV1,PV4,PV6,     
FCA0,FCA1,FCA4,FCA6,FCAAv,Zarr_Adj,DissConst,kCa,R,n,COAlk,B,TerrBioCPgC,Cstock1,Cstock2,
Raupach,k1A,k2A,WCARBs,BSILs,WSILs,volcs1,RiverP_mols,weathd13C,volcd13C,TerrRC,TerrBioSC,
TerrBioRC,BioSCF,BioRCF,Parr,Alkarr,Carr,SCarr,Fearr,Oarr,Siarr,RCarr,AtCO2,SCAt,RCAt,AtO2,
Sal0,Sal1,Sal2,Sal3,Sal4,Sal5,Sal6,TC0,TC1,TC2,TC3,TC4,TC5,TC6,rcp,rcpFe,rcpO,rcpSi,rphos,
SedCstock,SedAlkstock,ContCstock,AtCO2baseHol,AnthStock,DeforestCPgC,AnthSC,AnthRC,FASR,FSAR,
FK,FKR,RCS1At,RCS_atom,RCSAdjust,Avogrado,RCS_bombs,RCD1At,RCD1,SA_Earth_cm2,Sarr,MolC,Peta,
Tera,Giga,Pbar,secsday,secsyr,PerC,Varr,Kelv,SWD,Varrat,CVA,CVC,Sstand,Rstand,CDepth,FDepth,
OrgDepth,SedPbar,H,A1_C,A2_C,A3_C,B1_C,B2_C,B3_C,A1_O,A2_O,A3_O,B1_O,B2_O,B3_O)
from DataWork import (AnthAtd13C_time,AnthAtd13CDat,AnthAtD14C_time,AnthAtD14CDat,
AnthEmitsDat_time,AnthEmitsDat,AnthAtCO2_time,AnthAtCO2Dat_AtCO2,
rcp_time,rcp26_FFems,rcp45_FFems,rcp60_FFems,rcp85_FFems,rcp26_LUCems,rcp45_LUCems,rcp60_LUCems,rcp85_LUCems,
rcp26_AtCO2,rcp45_AtCO2,rcp60_AtCO2,rcp85_AtCO2,rcp26_SST,rcp45_SST,rcp60_SST,rcp85_SST,
GLODAP_time,GLODAP_DIC,GLODAP_Alk,GLODAP_Phos,GLODAP_d13C,GLODAP_D14C,GLODAP_CO23,
SD_GLODAP_DIC,SD_GLODAP_Alk,SD_GLODAP_Phos,SD_GLODAP_d13C,SD_GLODAP_D14C,ssp85_LUCems,ssp85_FFems,co2_add,CESM,Psi2_CESM)

import warnings
warnings.simplefilter("error")



##--------CONTROLS-------------------------------------------------------------
## User adjusted options, use these to edit model settings

Charting='off'# switches charting on or off (turn off for batch runs)
BatchRun='off' #turn off for manual model runs (ie those runs operated from this module)
StoreResults='off' # turn on to store results for starting data for future runs (turn off for batch runs) 
feedback='off'
j=0
y_add=0
x_add=0

TerrestrialGeo=1 #1 is on, 0 is off - switch on volcanic emissions and weathering sink
TerrestrialBios=1 #1 is on, 0 is off - switch on terrestrial biosphere
Rivers=1 #1 is on, 0 is off - switch on river fluxes of carbon, alkalinity and phosphorus
Carbonate=1 #1 is on, 0 is off - switch on carbonate chemistry
AnthEmits=1 #1 is on, 0 is off - switch on anthropogenic emissions
BombRadiocarbon=0 #1 is on, 0 is off - switch on bomb radiocarbon
## Remember - urn off 'AnthEmits' and 'BombRadiocarbon' when doing paleo runs
# with a time horizon that is not 'Anthrop' i.e. user-specified in years


# RCP settings, for modern scenarios
rcpname='rcp60' # Also need to adjust FFems, DeforestC, SST scenarios immediately below
rcp_FFems=ssp85_FFems
rcp_DeforestC=ssp85_LUCems
rcpSST=rcp85_SST
co2_add=co2_add
CESM=CESM
Psi2_CESM=Psi2_CESM

# rcp names ='rcp26','rcp45','rcp60','rcp85'

#Time calculations for modern day scenarios
AnthEmitStart=1751
AnthEmitEnd=2100
Anthrop=AnthEmitEnd-AnthEmitStart
BombStart=1954-AnthEmitStart
BombEnd=1963-AnthEmitStart

# Model run time (years) - important setting
interval = 1#Anthrop-85 # Model run time in years. Enter a time period in years,
# or if running the modern day emissions scenario, enter "Anthrop"

# Filenames for batch runs, scenarios and sensitivity tests - don't need to edit here
nh=0 # file counter, for naming of batch files
expname='n'# Setup for batch runs

##-----------------------------------------------------------------------------
## Generally do not need to vary settings below this line


## Other setup
dt =1 # model time step
TS=dt*secsyr # Time step calculation, converts seconds into years

# Establish the model----------------------------------------------------------
def SCPM_run(nh,expname,Psi1,Psi2,gamma1,gamma2,Zed,bScale,FCA0,FCA1,FCA4,FCA6,
             PV0,PV1,PV4,PV6,kCa,n,TerrBioCPgC,TC0,TC1,TC2,TC3,TC4,TC5,TC6,Sal0,Sal1,
             Sal4,Sal6,H,Varr,Parr,Siarr,Alkarr,SCarr,RCarr,Fearr,Oarr,Carr,AtCO2,SCAt, 
             RCAt,Cstock1,Cstock2,SedCstock,ContCstock,SedAlkstock,AnthStock,volcs1,weathd13C,AtO2):
     
    #Set up vectors
    time = []
    Psi1vec=[]
    Psi2vec=[]
    gamma1vec = []
    gamma2vec = []
    Pvec = []
    Cvec = []
    Alkvec = []
    pCO2vec = []
    pCO2Ovec=[]
    AtCO2vec=[]    
    SCvec = []
    SCAtvec=[]    
    RCvec = []
    RCAtvec=[]    
    Ovec = []
    Fevec = []
    Sivec = []
    AtO2vec=[]
    H2CO3vec=[]
    HCO3vec=[]
    CO23vec=[]
    DICvec=[]    
    AnthEmitvec=[]
    AnthStockvec=[]
    DeforestCvec=[]
    TCAtvec=[]
    TCTerrvec=[]
    TCSedvec=[]
    TotCarb_PgCvec=[]
    Cstock_PgCvec=[]
    Cstock_netvec=[]
    Cstock1_PgCvec=[]
    Cstock2_PgCvec=[]
    SedCstock_PgCvec=[]
    ContCstock_PgCvec=[]    
    TotCarbOc_PgCvec=[]
    CarbOc_PgCvec=[]
    TotCarbAt_PgCvec=[]
    OcAlkvec=[]
    SedAlkvec=[]
    TAlkvec=[]
    OCfluxvec=[]
    SumOCfluxvec=[]
    AtCfluxvec=[] 
    OCvec=[]
    Oargvec=[]
    pHvec=[]
    Hvec=[]
    K0vec=[]
    K1vec=[]
    LNK0vec=[]
    LNK1vec=[]
    LNK2vec=[]
    LNKspvec=[]
    DissCvec=[]
    NetSedCvec=[]
    SedCaCO3vec=[]
    results=[]
    SedCfluxvec=[]
    SumSedCfluxvec=[]
    SedDissCfluxvec=[]
    BioCfluxvec=[]
    BioC_Outfluxvec=[]
    Caflux_Outvec=[]
    DissCfluxvec=[]
    T1Cfluxvec=[]
    T2Cfluxvec=[]
    E1Cfluxvec=[]
    E2Cfluxvec=[]
    RivCfluxvec=[]
    VolcCfluxvec=[]
    WeathCfluxvec=[]
    TerrBioCfluxvec=[]
    CFertfluxvec=[]
    Respirefluxvec=[]
    RiverAlkvec=[]
    SedAlkfluxvec=[]
    RivCarbfluxvec=[]
    
    # Setup time series loop
    for i in range(interval):
        print(AtCO2)
        if j==0:
            AtCO2=4.146438732277601957e-04
        init_co2=4.146438732277601957e-04
    ## Human emissions and bomb radiocarbon activation (if switched on)--------
        AtCO2=AtCO2+(co2_add[j]+y_add+x_add)*1e-6
         #Industrial emissions - sources historical data and RCP projections
        AnthEmit1=0 # baseline setting
        if AnthEmits==1:
            AnthEmit1=rcp_FFems[-85+i+j]
        AnthEmit=(AnthEmit1*Peta/MolC)/secsyr # Industrial emissions
        
        # Deforestation emissions - sources historical data and RCP projections
        DeforestC=0 # baseline setting
        if AnthEmits==1:
            DeforestC=rcp_DeforestC[-85+i+j]
        DeforestC=DeforestC*Peta/MolC/secsyr # Deforestation emissions
        
        AnthStart=AnthEmitStart
        AnthStock=(AnthStock)+dt*secsyr*(AnthEmit+DeforestC) # Cumulative anthropogenic emissions
        
        # Bomb radiocarbon activation
        Bomb14C=RCS_bombs if i in range(BombStart,BombEnd) else 0
        
        #Anthropogenic emissions and carbon isotopes
        AnthSC1=AnthEmit*AnthSC
        AnthRC1=AnthEmit*AnthRC
        DeforestSC=(DeforestC*TerrBioSC)
        DeforestRC=(DeforestC*TerrBioRC)

    ## Implement matrices and other inputs-------------------------------------

        # For batch runs , parameters need to be implemented here
        
        # Temperature sensitivity
        TCadjust=0
        if AnthEmits==1:
            TCadjust=rcpSST[-85+i+j]
        
        # Temperature
        TC = np.zeros([7,1])
        TC[0,0]=TC0+TCadjust # RCP SST perturbations (warming) or LGM cooling
        TC[1,0]=TC1+TCadjust # RCP SST perturbations (warming) or LGM cooling
        TC[2,0]=TC2
        TC[3,0]=TC3
        TC[4,0]=TC4
        TC[5,0]=TC5
        TC[6,0]=TC6+TCadjust # RCP SST perturbations (warming) or LGM cooling
        
        # Salinity
        Sal=np.zeros([7,1]) 
        Sal[0,0]=Sal0
        Sal[1,0]=Sal1
        Sal[2,0]=Sal2
        Sal[3,0]=Sal3
        Sal[4,0]=Sal4
        Sal[5,0]=Sal5
        Sal[6,0]=Sal6
        
        # Rain ratio
        FCA=np.ones([7,1]) 
        FCA[0,0]=FCA0
        FCA[1,0]=FCA1
        FCA[4,0]=FCA4
        FCA[6,0]=FCA6

        # Convert degrees C to kelvin
        TK = (TC+Kelv) 
        
        # Piston velocity
        FXarr=makeFXarr(PV0,PV1,PV4,PV6,secsday)*Sarr
  
        # Multiply ocean physical parameters by their respective matrix
        T1mat = makeT1mat(Psi1,Varr,Fract)*Psi1/Varr
        T2mat = makeT2mat(Psi2,Varr)*Psi2/Varr
        E1mat=makeE1mat(gamma1,Varr)*gamma1/Varr
        E2mat = makeE2mat(gamma2,Varr)*gamma2/Varr
        PhysMat=T1mat+T2mat+E1mat+E2mat
        
    ## Biological pump---------------------------------------------------------
        
        # Implement scalar on Zed biological flux parameter
        Zarr=makeZarr(Zed)*Zarr_Adj 
        OrgCarr=Zarr*Sarr/secsyr
        OrgSCarr=OrgCarr*BioSCF*Sstand
        OrgRCarr=OrgCarr*BioRCF*Rstand
        OrgParr=OrgCarr*rphos
        
        # Implement depth scalar function for biological flux (Martin et al, 1987)
        CDarr=(CDepth/OrgDepth)**-bScale 
        FDarr=(FDepth/OrgDepth)**-bScale
        CDarr=np.nan_to_num(CDarr)
        FDarr=np.nan_to_num(FDarr)                
        
        # Implement biological pump
        BioOut=makeBiomat_out(Zed, Varr)*FDarr
        BioIn=makeBiomat_in(Zed, Varr)*CDarr
        BioP = np.dot(BioIn+BioOut,OrgParr)/Varr 
        BioC = np.dot(BioIn+BioOut,OrgCarr)/Varr
        BioFe=np.dot(BioP,rcpFe)
        BioSi=np.dot(BioP,rcpSi)
        BioO=np.dot(BioP,rcpO)
        
        # Carbon isotopes
        BioSC = np.dot(BioIn+BioOut,OrgSCarr)/Varr
        BioRC = np.dot(BioIn+BioOut,OrgRCarr)/Varr

    ## Carbonate pump----------------------------------------------------------
        
    # Rain ratio x organic carbon flux at 100m depth
        CaCOflux=FCA*Zarr       
        CaCO3flux=(CaCOflux*Sarr)/Varr/secsyr   
        Caflux_out1=-CaCO3flux
        
    ## pCO2 and carbonate calcs------------------------------------------------
    
        # The K coefficients, references shown
        IonS = 19.924*Sal/(1000-1.005*Sal) # Salinity ion
        BT = 1.179e-5*Sal # Total Boron
        Ca=0.01028*Sal/35 # Calcium
        K0 = np.exp((A1_C + A2_C * (100.0/((TC+Kelv))) + A3_C*np.log((TC+Kelv)/100.0) + (Sal) * (B1_C + B2_C *((TC+Kelv)/100.0) + B3_C * (((TC+Kelv)/100.0)**2)))); # Weiss (1974) in mols kg-1 atm-1
        K1 = 10**-((3633.86/(TC+Kelv)) - 61.2172 + (9.67770 * np.log(TC+Kelv)) - (0.011555 * Sal) + (0.0001152 * (Sal**2))); #Lueker et al (2000) in mol kg-1
        K2 = 10**-((471.78/((TC+Kelv))) + 25.9290 - (3.16967 * np.log((TC+Kelv))) - (0.01781 * Sal) + (0.0001122 * (Sal**2))); # Lueker et al (2000) in mol kg-1
        KB = np.exp((-8966.9-2890.53*(Sal**0.5)-77.942*Sal+1.728*(Sal**1.5)-0.0996*(Sal**2))/(TC+Kelv)+(148.0248+137.1942*(Sal**0.5)+(1.62142*Sal))+(-24.4344-25.085*(Sal**0.5)-0.2474*Sal)*np.log((TC+Kelv))+0.053105*(Sal**0.5)*(TC+Kelv)) #Dickson 1990
        KW = np.exp(148.96502+(-13847.26/(TC+Kelv))-(23.6521*np.log((TC+Kelv)))+((Sal**0.5)*(-5.977+118.67/(TC+Kelv)+1.0495*np.log((TC+Kelv))))-0.01615*Sal) #Millero (1995)
        KS = np.exp(-8904.2/(TC+Kelv)+117.4-19.334*np.log((TC+Kelv))+(-458.79/(TC+Kelv)+3.5913)*IonS**0.5+(188.74/(TC+Kelv)-1.5998)*IonS+(-12.1652/(TC+Kelv)+0.07871)*IonS*IonS)*(1-0.001005*Sal)
        KP1 = np.exp(-4576.752/(TC+Kelv)+115.54-18.453*np.log((TC+Kelv))+(-106.736/(TC+Kelv)+0.69171)*Sal**0.5+(-0.65643/(TC+Kelv)-0.01844)*Sal) #Yao and Millero (1995)
        KP2 = np.exp(-8814.715/(TC+Kelv)+172.1033-27.927*np.log((TC+Kelv))+(-160.34/(TC+Kelv)+1.3566)*Sal**0.5+(0.37335/(TC+Kelv)-0.05778)*Sal) #Yao and Millero (1995)
        KP3 = np.exp(-3070.75/(TC+Kelv)-18.126+(17.27039/(TC+Kelv)+2.81197)*Sal**0.5+(-44.99486/(TC+Kelv)-0.09984)*Sal) #Yao and Millero (1995)
        KF = np.exp(1590.2/(TC+Kelv)-12.641+1.525*IonS**0.5)*(1-0.001005*Sal)
        Ksp1 = np.exp(-395.8293+(6537.773/(TC+Kelv))+71.595*np.log((TC+Kelv))-0.17959*(TC+Kelv)+(-1.78938+410.64/(TC+Kelv)+0.0065453*(TC+Kelv))*Sal**0.5-0.17755*Sal+0.0094979*Sal**(3.0/2.0)) # Mucci (1983)
        Karg1 = np.exp(-395.9180+(6685.079/(TC+Kelv))+71.595*np.log((TC+Kelv))-0.17959*(TC+Kelv)+(-0.157481+202.938/(TC+Kelv)+0.003978*(TC+Kelv))*Sal**0.5-0.23067*Sal+0.0136808*Sal**(3.0/2.0)) # Mucci (1983)
        K02= np.exp((A1_O + A2_O * (100.0/((TC+Kelv))) + A3_O*np.log((TC+Kelv)/100.0) + (Sal) * (B1_O + B2_O *((TC+Kelv)/100.0) + B3_O * (((TC+Kelv)/100.0)**2)))); # Weiss (1974) in mols kg-1 atm-1
     
        # LNK's for checking against literature values
        LNK0 =np.log(K0)/np.log(0.1)
        LNK1=np.log(K1)/np.log(0.1)
        LNK2=np.log(K2)/np.log(0.1)
        LNKsp=np.log(Ksp1)/np.log(0.1) 
        LNKarg=np.log(Karg1)/np.log(0.1) 
        
        # Pressure corrections for Ksp as per Millero (1983)
        dVC=-65.28+0.397*TC-0.005155*TC**2+(19.816-0.0441*TC-0.00017*TC**2)*(Sal/35)**0.5
        dKC=0.01847+0.0001956*TC-0.000002212*TC**2+(-0.03217-0.0000711*TC-0.000002212*TC**2)*(Sal/35)**0.5
        dVA=-65.50+0.397*TC-0.005155*TC**2+(19.82-0.0441*TC-0.00017*TC**2)*(Sal/35)**0.5
        dKA=0.01847+0.0001956*TC-0.000002212*TC**2+(-0.03217-0.0000711*TC-0.000002212*TC**2)*(Sal/35)**0.5
        Kspres = (-dVC+0.5*dKC*Pbar)*Pbar/(R*(TC+Kelv))
        Ksp=Ksp1*np.exp(Kspres)
        Kargpres= (-dVA+0.5*dKA*Pbar)*Pbar/(R*(TC+Kelv))
        Karg=Karg1*np.exp(Kargpres)
        
        # For sediments, at the bottom of the the abyssal box (box 5)
        KspresSed = (-dVC[5,0]+0.5*dKC[5,0]*SedPbar)*SedPbar/(R*(TC[5,0]+Kelv))
        KspSed=Ksp1[5,0]*np.exp(KspresSed)        
        
        # Create pH components (Follows et al, 2006)
        BAlk=BT*KB/(H+KB)
        denom=H*H*H+(KP1*H*H)+(KP1*KP2*H)+(KP1*KP2*KP3)
        SiAlk=(Siarr/SWD)*KS/(KS+H)
        H3PO4g=((Parr/SWD)*H*H*H)/denom
        H2PO4g=((Parr/SWD)*KP1*H*H)/denom
        HPO4g=((Parr/SWD)*KP1*KP2*H)/denom
        PO4g=((Parr/SWD)*KP1*KP2*KP3)/denom
        fg = -BAlk-(KW/H)+H-HPO4g-2*PO4g+H3PO4g-SiAlk
        
        # Calculate pCO2 and carbon chemistry (Follows et al, 2006)
        CAlk=Alkarr/SWD+fg #Alkarr converted back into mol kg-1 via sea water density (SWD)
        gamma=Carr/SWD/CAlk        
        DM=(1-gamma)*(1-gamma)*K1*K1-4*K1*K2*(1-2*gamma)
        H=0.5*((gamma-1)*K1+(DM)**0.5)
        pCO2 = (Carr/SWD/K0)*((H**2)/((H**2)+K1*H+K1*K2)) # Carr converted back into mol kg-1 via SWD
        pH = np.log(H)/np.log(0.1)
        H2CO3 = pCO2*K0 # mol kg-1 
        HCO3 = (H2CO3*K1)/H # mol kg-1
        CO23 = (HCO3*K2)/H # mol kg-1
        DIC = H2CO3+HCO3+CO23 # mol kg-1
        
        BAlk=BT*KB/(H+KB)
        denom=H*H*H+(KP1*H*H)+(KP1*KP2*H)+(KP1*KP2*KP3)
        SiAlk=(Siarr/SWD)*KS/(KS+H)
        H3PO4g=((Parr/SWD)*H*H*H)/denom
        H2PO4g=((Parr/SWD)*KP1*H*H)/denom
        HPO4g=((Parr/SWD)*KP1*KP2*H)/denom
        PO4g=((Parr/SWD)*KP1*KP2*KP3)/denom
        fg = -BAlk-(KW/H)+H-HPO4g-2*PO4g+H3PO4g-SiAlk
        
        # Calculate pCO2 and carbon chemistry (Follows et al, 2006)
        CAlk=Alkarr/SWD+fg #Alkarr converted back into mol kg-1 via sea water density (SWD)
        gamma=Carr/SWD/CAlk        
        DM=(1-gamma)*(1-gamma)*K1*K1-4*K1*K2*(1-2*gamma)
        H=0.5*((gamma-1)*K1+(DM)**0.5)
        pCO2 = (Carr/SWD/K0)*((H**2)/((H**2)+K1*H+K1*K2)) # Carr converted back into mol kg-1 via SWD
        pH = np.log(H)/np.log(0.1)
        H2CO3 = pCO2*K0 # mol kg-1 
        HCO3 = (H2CO3*K1)/H # mol kg-1
        CO23 = (HCO3*K2)/H # mol kg-1
        DIC = H2CO3+HCO3+CO23 # mol kg-1
        
    # Carbonate dissolution calculations---------------------------------------

        # Omegas for dissolution calculations
        Ca=0.01028*Sal/35.0
        OmegaC=np.minimum(CO23*Ca/Ksp,1)
        OmegaArg=np.minimum(CO23*Ca/Karg,1)
        OmegaCSed=np.minimum(CO23[5,0]*Ca/KspSed,1)

        # Water column
        DissC=np.zeros([7,1])
        DissC = (CO23*Ca)*(SWD*1e3)*kCa*(1-OmegaC)**n
        DissC=(DissC*1e-3*PerC+DissConst)*Carbonate # converted from mmol/m3 back to mol/m3 and then adjusted for % C in CaCO3
        
        # Sediments
        SedDiss=np.zeros([7,1])
        SedDiss[5,0]=1 # to isolate sediment dissolution flux to the abyssal box
        SedDissC=(CO23*Ca)*(SWD*1e3)*kCa*(1-OmegaCSed)**n
        SedDissC1=(SedDissC*1e-3*PerC+DissConst)*SedDiss
        
        # Total carbonate flux (sinking particles, dissolution)
        NetCflux=(Caflux_out1+DissC+SedDissC1)*Carbonate
        NetAlkflux=(Caflux_out1*2.0+DissC*2.0+SedDissC1*2.0)*Carbonate        
            
    ## Air-sea gaseous fluxes--------------------------------------------------

        # CO2 fluxes
        cflux1=((SWD*K0*FXarr*(AtCO2-pCO2)))# SWD converts mol/kg into mol/m3
        cflux=cflux1/Varr
        Atcflux=-sum(cflux1)/Varrat   

        # d13C air-sea fractionation factors
        FSA = np.zeros([7,1])
        FSA=(-9.866/TK+1.02412) # Mook (1974)              
        FAS = np.zeros([7,1])
        FAS=(-0.373/TK+1.00019) # Mook (1974)        
    
        # Air-sea flux of 13C
        SCPCO2 = (SWD*K0*FXarr*FK*((FAS*(SCAt/AtCO2))*AtCO2-(FSA*(SCarr/Carr))*pCO2))
        Scflux=SCPCO2/Varr # Ocean boxes
        AtSCflux=-sum(SCPCO2)/Varrat # Atmosphere
        
        # Air-sea flux of 14C
        RCPCO2 =SWD*K0*FXarr*FKR*((FASR*(RCAt/AtCO2))*AtCO2-(FSAR*(RCarr/Carr)*pCO2))
        Rcflux=RCPCO2/Varr # Ocean boxes   
        AtRCflux=-sum(RCPCO2)/Varrat # Atmosphere
        
        # Air-sea oxygen flux
        OO2=SWD*K02*AtO2
        Oflux1=(FXarr*(OO2-Oarr))# SWD converts mol/kg into mol/m3
        Oflux=Oflux1/Varr
        AtOflux=-sum(Oflux1)/Varrat  

    ## Atmospheric source of radiocarbon
        RCS_mol=RCS_atom/Avogrado
        RCS1At =RCS_mol*SA_Earth_cm2

    ## Weathering, river fluxes and volcanic emissions-------------------------
        
        # As per Toggweiler (2007) only silicate weathering is a sink of CO2 from the atmosphere
        # Weathering of carbonate rocks is a source of carbon to the low latitude surface ocean via rivers
        RVCARB=WCARBs*AtCO2
        RVSIL=(BSILs+WSILs*AtCO2)
        weaths=RVSIL*Varr[0,0] #silicate rock weathering sink of CO2, carbonate weathering is a source of carbon in rivers
        
        # Volcanic carbon emissions
        volcs=weaths # volcanic emissions in step with silicate weathering as per Toggweiler (2007)
        # unless anthropocene scenario, hardwired estimate 
        if AnthEmits==1:
            volcs=volcs1
        else: 
            volcs=weaths
        
        # net source/sink of terrestrial carbon
        TerrC=(volcs-weaths)/Varrat 
        TCflux=TerrC*TerrestrialGeo
        
        # weathering and volcanism 13C and 14C fluxes
        weathd13C=weathd13C
        TerrSC=(-weaths*weathd13C+volcs*volcd13C)/Varrat     
        TSCflux=TerrSC*TerrestrialGeo
        TRCflux=TerrC*TerrRC*TerrestrialGeo

        # River fluxes
        RiverCflux=np.zeros([7,1])
        RiverCflux[0,0]=(RVCARB+RVSIL) # mol/m3 Incorporates source of carbonate weathering        
        RiverAlkflux=RiverCflux*2.0 #Alk:C ratio 2:1 as per Toggweiler (2007)
        RiverPflux=np.zeros([7,1])
        RiverPflux[0,0]=RiverP_mols/Varr[0,0]
        PSedflux=np.zeros([7,1])        
        PSedflux[5,0]=RiverP_mols/Varr[5,0]

    ## Terrestrial biosphere---------------------------------------------------
        
        # Carbon fertilisation and respiration
        Ratio=np.maximum(AtCO2/AtCO2baseHol,0) # AtCO2 versus a baseline level
        TerrBioC=TerrBioCPgC*Peta/MolC/secsyr
        CFert=(TerrBioC*(1+B*np.log(Ratio)))*TerrestrialBios 
        Respire1=k1A*Cstock1*TerrestrialBios        
        Respire2=k2A*Cstock2*TerrestrialBios 
        Respire=Respire1+Respire2
        
        # Stocks
        Cstock1=Cstock1+dt*secsyr*(CFert*Raupach-Respire1)*TerrestrialBios
        Cstock2=Cstock2+(dt*secsyr*(CFert*(1-Raupach)-Respire2)-(dt*secsyr*DeforestC*AnthEmits))*TerrestrialBios        
        Cstock=Cstock1+Cstock2

    ## Step forward model calculations-----------------------------------------
      
        # Model equations
        
        if feedback='on':    
            BioC[1,0]=(np.log(AtCO2*1e6)*0.01654927-0.13269258)/(np.log(init_co2*1)*0.01654927-0.13269258)*-2.87085258e-10 #BDIC[0]
            BioA=np.zeros((7,1))
            BioA=BioC*-16/106
            BioA[1,0]=(np.log(AtCO2*1e6)*-0.0061603+0.04017557)/np.abs(np.log(init_co2)*-0.0061603+0.04017557)*2.87085258e-10*16/106 #BALK[0]
        else:
            BioA=np.zeros((7,1))
            BioA=-16/106*BioC #BALK[0]
        
        # Ocean boxes
        Parr = Parr +dt*secsyr*(np.dot(PhysMat, Parr)+BioP+RiverPflux*Rivers-PSedflux*Rivers)
        Carr = Carr + dt*secsyr*(np.dot(PhysMat, Carr)+BioC+cflux+NetCflux+RiverCflux*Rivers)
        Alkarr = Alkarr + dt*secsyr*(np.dot(PhysMat, Alkarr)+NetAlkflux+BioA+RiverAlkflux*Rivers)
        Fearr = Fearr + dt*secsyr*(np.dot(PhysMat, Fearr)+BioFe)
        Siarr = Siarr + dt*secsyr*(np.dot(PhysMat, Siarr)+BioSi)
        Oarr = Oarr + dt*secsyr*(np.dot(PhysMat, Oarr)+BioO)
        SCarr = SCarr+ dt*secsyr*(np.dot(PhysMat, SCarr)+BioSC+Scflux+NetCflux*Sstand+RiverCflux*Sstand*Rivers)
        SCratio=SCarr/Carr
        RCarr = RCarr + dt*secsyr*(np.dot(PhysMat, RCarr)+BioRC+Rcflux+NetCflux*(RCarr/Carr)-(RCD1*RCarr*Varr)/Varr)
        
        AC=(np.dot(PhysMat, Carr))
        AA=(np.dot(PhysMat, Alkarr))
  
        # Atmosphere
        AtCO2 = AtCO2 + dt*secsyr*(Atcflux+TCflux-((CFert-Respire)/Varrat)+((AnthEmit+DeforestC)/Varrat)*AnthEmits)
        pCO2a=np.append(pCO2,AtCO2) # create an array of all pCO2
        SCAt = SCAt+dt*secsyr*(AtSCflux+TSCflux-((CFert-Respire)/Varrat)*TerrBioSC*TerrestrialBios+(AnthSC1+DeforestSC)/Varrat*AnthEmits)
        RCAt = RCAt+dt*secsyr*(AtRCflux+RCS1At/Varrat-RCD1At*RCAt*Varrat/Varrat-TRCflux-((CFert-Respire)/Varrat)*TerrBioRC*TerrestrialBios+((AnthRC1+DeforestRC)/Varrat)*AnthEmits+(Bomb14C/Varrat)*BombRadiocarbon)
        AtO2=AtO2
 
    ## Carbon stocks in PgC----------------------------------------------------

        # Total carbon in each reservoir
        CarbOc_PgC=(Carr*Varr*MolC)/Peta # Each individual ocean box
        TotCarbOc_PgC=(sum(Carr*Varr)*MolC)/Peta # Sum of ocean boxes
        TotCarbAt_PgC=(AtCO2*Varrat*MolC)/Peta # Atmosphere
        SedCstock=SedCstock-sum(NetCflux*secsyr*Varr)-sum(BioC*secsyr*Varr) # Marine sediments
        SedCstock_PgC=SedCstock*MolC/Peta
        ContCstock=ContCstock-(dt*secsyr*volcs)-(dt*secsyr*RVCARB*Varr[0,0]) #Continents
        ContCstock_PgC=(ContCstock*MolC)/Peta
        
        # Terrestrial biosphere
        Cstock_PgC=(Cstock*MolC)/Peta
        Cstock1_PgC=(Cstock1*MolC)/Peta
        Cstock2_PgC=(Cstock2*MolC)/Peta

        # Total carbon in earth system
        TotCarb_PgC=TotCarbOc_PgC+TotCarbAt_PgC+Cstock_PgC+SedCstock_PgC+ContCstock_PgC   
        
        # Alkalinity stocks
        OcAlk_Tmol=sum(Alkarr*Varr)/Tera
        SedAlkstock=SedAlkstock-sum(NetAlkflux*secsyr*Varr)
        SedAlkstock_Tmol=SedAlkstock/Tera
        TotAlk=OcAlk_Tmol+SedAlkstock_Tmol
        
    ##Carbon fluxes in PgC/yr--------------------------------------------------
        
        SedCflux_PgC=(secsyr*NetCflux*Varr*MolC)/Peta
        SumSedCflux_PgC=sum(SedCflux_PgC)  
        BioCflux_PgC=(secsyr*BioC*Varr*MolC)/Peta
        Caflux_OutPgC=(secsyr*Caflux_out1*Varr*MolC)/Peta # Carbonate fluxes from the surface boxes
        DissC_PgC=(secsyr*DissC*Varr*MolC)/Peta        
        SedDissC_PgC=(secsyr*SedDissC1*Varr*MolC)/Peta
        SumSedDissC_PgC=sum(SedDissC_PgC)# sum out zero values 
        T1Cflux_PgC=(secsyr*np.dot(T1mat, Carr)*Varr*MolC)/Peta
        T2Cflux_PgC=(secsyr*np.dot(T2mat, Carr)*Varr*MolC)/Peta
        E1Cflux_PgC=(secsyr*np.dot(E1mat, Carr)*Varr*MolC)/Peta
        E2Cflux_PgC=(secsyr*np.dot(E2mat, Carr)*Varr*MolC)/Peta    
        VolcCflux_PgC=(secsyr*volcs*MolC)/Peta*TerrestrialGeo
        WeathCflux_PgC=(secsyr*-weaths*MolC)/Peta*TerrestrialGeo
        RivCflux=sum(RiverCflux*Varr)
        RivCflux_PgC=((secsyr*RivCflux*MolC)/Peta)*Rivers
        RiverCarbflux_PgC=(secsyr*RVCARB*Varr[0,0]*MolC)/Peta*TerrestrialGeo # carbonate weathering influx
        CFertflux_PgC=(-CFert*MolC*dt*secsyr)/Peta
        Respireflux_PgC=(Respire*MolC*dt*secsyr)/Peta        
        TerrBioCflux_PgC=(CFertflux_PgC+Respireflux_PgC)# baseline input, not used       
        
        # Air-sea fluxes of carbon - note multiplied back out by volumes as the model diagnostics are per m3
        SumOCflux=(sum(cflux1)*dt*secsyr*MolC)/Peta
        OCflux=(cflux1*dt*secsyr*MolC)/Peta
        AtCflux=(Atcflux*Varrat*dt*secsyr*MolC)/Peta

        # Fluxes of alkalinity
        RiverAlk_Tmol=((secsyr*RiverAlkflux[0,0]*Varr[0,0])/Tera)*Rivers
        SedAlkflux=sum((NetAlkflux*secsyr*Varr)/Tera)
    
    ## Results processing and output-------------------------------------------
    
    ## Store output in a time series
        if np.mod(i,1)==0:
            time = np.append(time,i)
            Psi1vec= np.append(Psi1vec, Psi1/1e6)
            Psi2vec=np.append(Psi2vec, Psi2/1e6)
            AnthEmitvec=np.append(AnthEmitvec,AnthEmit*secsyr*MolC/Peta)
            AnthStockvec=np.append(AnthStockvec,AnthStock*MolC/Tera)
            DeforestCvec=np.append(DeforestCvec,DeforestC*secsyr*MolC/Peta)
            Pvec = np.append(Pvec,Parr)
            Cvec = np.append(Cvec, Carr)
            AtCO2vec= np.append(AtCO2vec, AtCO2)
            AtO2vec= np.append(AtO2vec, AtO2)
            Alkvec = np.append(Alkvec, Alkarr)        
            SCvec = np.append(SCvec,SCarr)
            SCAtvec=np.append(SCAtvec,SCAt)
            RCvec = np.append(RCvec,RCarr)    
            RCAtvec = np.append(RCAtvec,RCAt)
            Fevec = np.append(Fevec, Fearr)
            Sivec = np.append(Sivec, Siarr)        
            Ovec = np.append(Ovec, Oarr)          
            Cstock_PgCvec=np.append(Cstock_PgCvec, Cstock_PgC)           
            Cstock1_PgCvec=np.append(Cstock1_PgCvec, Cstock1_PgC)
            Cstock2_PgCvec=np.append(Cstock2_PgCvec, Cstock2_PgC)            
            SedCstock_PgCvec=np.append(SedCstock_PgCvec, SedCstock_PgC)                        
            ContCstock_PgCvec=np.append(ContCstock_PgCvec, ContCstock_PgC)
            DICvec=np.append(DICvec,DIC)
            TotCarb_PgCvec=np.append(TotCarb_PgCvec, TotCarb_PgC)
            TotCarbOc_PgCvec=np.append(TotCarbOc_PgCvec, TotCarbOc_PgC)
            CarbOc_PgCvec=np.append(CarbOc_PgCvec, CarbOc_PgC)
            TotCarbAt_PgCvec=np.append(TotCarbAt_PgCvec, TotCarbAt_PgC)          
            TAlkvec=np.append(TAlkvec, TotAlk)
            OcAlkvec=np.append(OcAlkvec, OcAlk_Tmol)
            SedAlkvec=np.append(SedAlkvec, SedAlkstock_Tmol)                           
            pCO2Ovec = np.append(pCO2Ovec, pCO2*1e6)
            pCO2vec = np.append(pCO2vec, pCO2a*1e6) # ppm for presentation purposes          
            H2CO3vec=np.append(H2CO3vec,H2CO3*1e6) # note this is umol/kg for presentation and comparison with data
            HCO3vec=np.append(HCO3vec,HCO3**1e6) # note this is umol/kg for presentation and comparison with data
            CO23vec=np.append(CO23vec,CO23*1e6) # note this is umol/kg for presentation and comparison with data
            OCvec=np.append(OCvec,OmegaC)
            Oargvec=np.append(Oargvec,OmegaArg)
            pHvec=np.append(pHvec,pH)
            DissCvec=np.append(DissCvec,DissC_PgC)
            NetSedCvec=np.append(NetSedCvec,NetCflux)
            SumOCfluxvec=np.append(SumOCfluxvec, SumOCflux)
            OCfluxvec=np.append(OCfluxvec, OCflux)     
            AtCfluxvec=np.append(AtCfluxvec, AtCflux)               
            SedCfluxvec=np.append(SedCfluxvec, SedCflux_PgC)
            SumSedCfluxvec=np.append(SumSedCfluxvec, SumSedCflux_PgC)
            BioCfluxvec=np.append(BioCfluxvec, BioCflux_PgC)
            Caflux_Outvec=np.append(Caflux_Outvec, Caflux_OutPgC)
            DissCfluxvec=np.append(DissCfluxvec,DissC_PgC)
            SedDissCfluxvec=np.append(SedDissCfluxvec,SumSedDissC_PgC)
            T1Cfluxvec=np.append(T1Cfluxvec, T1Cflux_PgC)
            T2Cfluxvec=np.append(T2Cfluxvec, T2Cflux_PgC)
            E1Cfluxvec=np.append(E1Cfluxvec, E1Cflux_PgC)
            E2Cfluxvec=np.append(E2Cfluxvec, E2Cflux_PgC)            
            VolcCfluxvec=np.append(VolcCfluxvec, VolcCflux_PgC)
            WeathCfluxvec=np.append(WeathCfluxvec, WeathCflux_PgC)
            RivCfluxvec=np.append(RivCfluxvec, RivCflux_PgC)
            RivCarbfluxvec=np.append(RivCarbfluxvec, RiverCarbflux_PgC)            
            TerrBioCfluxvec=np.append(TerrBioCfluxvec, TerrBioCflux_PgC)
            CFertfluxvec=np.append(CFertfluxvec, CFertflux_PgC)
            Respirefluxvec=np.append(Respirefluxvec, Respireflux_PgC)                        
            RiverAlkvec=np.append(RiverAlkvec, RiverAlk_Tmol)
            SedAlkfluxvec=np.append(SedAlkfluxvec, SedAlkflux)
            
            
    # Prepare output for charting
    nt = time.shape[0]
    Phos = (np.reshape(Pvec,(nt,7)))
    Carb = (np.reshape(Cvec,(nt,7)))
    AtCO2=(np.reshape(AtCO2vec,(nt,1)))   
    AtO2=(np.reshape(AtO2vec,(nt,1)))  
    Alk = (np.reshape(Alkvec,(nt,7)))
    SCarb = (np.reshape(SCvec,(nt,7)))
    SCAt=(np.reshape(SCAtvec,(nt,1)))    
    RCarb = (np.reshape(RCvec,(nt,7)))
    RCAt=(np.reshape(RCAtvec,(nt,1)))    
    Fe = (np.reshape(Fevec,(nt,7)))
    Si = (np.reshape(Sivec,(nt,7)))
    O = (np.reshape(Ovec,(nt,7)))    
    Cstock_PgC=(np.reshape(Cstock_PgCvec,(nt,1)))    
    Cstock1_PgC=(np.reshape(Cstock1_PgCvec,(nt,1)))
    Cstock2_PgC=(np.reshape(Cstock2_PgCvec,(nt,1))) 
    TotCarbOc_PgC=(np.reshape(TotCarbOc_PgCvec,(nt,1)))
    CarbOc_PgC=(np.reshape(CarbOc_PgCvec,(nt,7)))
    TotCarbAt_PgC=(np.reshape(TotCarbAt_PgCvec,(nt,1)))
    SedCstock_PgC=(np.reshape(SedCstock_PgCvec,(nt,1)))
    ContCstock_PgC=(np.reshape(ContCstock_PgCvec,(nt,1)))
    TotCarb_PgC=(np.reshape(TotCarb_PgCvec,(nt,1)))
    AnthEmitOut=(np.reshape(AnthEmitvec,(nt,1))) # converted to PgC from TgC
    AnthStockOut=(np.reshape(AnthStockvec,(nt,1)))/1000 # converted to PgC from TgC
    DeforestCOut=(np.reshape(DeforestCvec,(nt,1)))
    pCO2a = (np.reshape(pCO2vec,(nt,8)))
    H2CO3out = (np.reshape(H2CO3vec,(nt,7)))
    HCO3out = (np.reshape(HCO3vec,(nt,7)))
    CO23out = (np.reshape(CO23vec,(nt,7)))
    DIC = (np.reshape(DICvec,(nt,7)))
    pH= (np.reshape(pHvec,(nt,7)))
    NetSedC=(np.reshape(NetSedCvec,(nt,7)))
    AtCflux=(np.reshape(AtCfluxvec,(nt,1)))
    OcCflux=(np.reshape(OCfluxvec,(nt,7)))
    SumOcCflux=(np.reshape(SumOCfluxvec,(nt,1)))    
    SedCflux=(np.reshape(SedCfluxvec,(nt,7)))
    SumSedCflux=(np.reshape(SumSedCfluxvec,(nt,1)))
    SedDissCflux=(np.reshape(SedDissCfluxvec,(nt,1)))
    BioCflux=(np.reshape(BioCfluxvec,(nt,7))) 
    Caflux_Out=(np.reshape(Caflux_Outvec,(nt,7)))
    DissCflux=(np.reshape(DissCvec,(nt,7)))
    T1Cflux=(np.reshape(T1Cfluxvec,(nt,7)))
    T2Cflux=(np.reshape(T2Cfluxvec,(nt,7)))
    E1Cflux=(np.reshape(E1Cfluxvec,(nt,7)))
    E2Cflux=(np.reshape(E2Cfluxvec,(nt,7)))
    VolcCflux=(np.reshape(VolcCfluxvec,(nt,1)))
    WeathCflux=(np.reshape(WeathCfluxvec,(nt,1)))
    RivCflux=(np.reshape(RivCfluxvec,(nt,1)))
    RivCarbCflux=(np.reshape(RivCarbfluxvec,(nt,1))) 
    TerrBioCflux=(np.reshape(TerrBioCfluxvec,(nt,1)))
    CFertflux=(np.reshape(CFertfluxvec,(nt,1)))
    Respireflux=(np.reshape(Respirefluxvec,(nt,1))) 
    RiverAlkflux_Tmol=(np.reshape(RiverAlkvec,(nt,1)))
    SedAlkflux_Tmol=(np.reshape(SedAlkfluxvec,(nt,1)))    
    TotalEmitOut=(AnthEmitOut+DeforestCOut)
  
    # Unit transformations for reporting/charting
    PhosOut=Phos*CVC
    CarbOut=Carb*CVC 
    AtCO2Out=AtCO2*1e6
    AtO2Out=AtO2
    AlkOut=Alk*CVC
    SCarbOut=((SCarb/Carb/Sstand)-1)*1000
    SCAtOut=((SCAt/AtCO2/Sstand)-1)*1000
    RCarbOut=((RCarb/Carb/Rstand)-1)*1000
    RCAtOut=((RCAt/AtCO2/Rstand)-1)*1000    
    FeOut=Fe*CVC
    SiOut=Si*CVC
    OOut=O*CVC
    DICOut=DICvec*1e6
  
    # Processed output for model-data file (equilibrium/end of simulation values)
    EqPhosOut=PhosOut[interval-1]-np.zeros([1,7])
    EqCarbOut=CarbOut[interval-1]-np.zeros([1,7])
    EqCarb=Carb[interval-1]-np.zeros([1,7]) # mol/m3 for CC diagram calcs
    EqAtCO2Out=AtCO2Out[interval-1]-np.zeros([1,1])
    EqAtO2Out=AtO2Out[interval-1]-np.zeros([1,1])
    EqAlkOut=AlkOut[interval-1]-np.zeros([1,7])
    Eqd13COut=SCarbOut[interval-1]-np.zeros([1,7])
    EqAtd13COut=SCAtOut[interval-1]-np.zeros([1,1])
    Eqd14COut=RCarbOut[interval-1]-np.zeros([1,7])
    EqAtd14COut=RCAtOut[interval-1]-np.zeros([1,1])
    EqSiOut=SiOut[interval-1]-np.zeros([1,7])
    EqFeOut=FeOut[interval-1]-np.zeros([1,7])
    EqOOut=OOut[interval-1]-np.zeros([1,7])
    EqCstock1_PgC=Cstock1_PgC[interval-1]-np.zeros([1,1])
    EqCstock2_PgC=Cstock2_PgC[interval-1]-np.zeros([1,1])
    EqCO23Out=CO23out[interval-1]-np.zeros([1,7])#LGMCO23
    RCAtOutav=np.average(RCAtOut)

    # Carbon reservoir stocks and fluxes at equilibrium/end of simulation
    # Many of these may not be needed, generally, but form the inputs
    # for a carbon cycle diagram (stocks and fluxes )
    EqTotCarbOc_PgC=TotCarbOc_PgC[interval-1]
    EqTotCarbAt_PgC=TotCarbAt_PgC[interval-1]
    EqSedCstock_PgC=SedCstock_PgC[interval-1]
    EqContCstock_PgC=ContCstock_PgC[interval-1]
    EqCstock_PgC=Cstock_PgC[interval-1]
    EqCarbOc_PgC=CarbOc_PgC[interval-1]
    EqTotCarb_PgC=TotCarb_PgC[interval-1] 
    EqAtCflux=AtCflux[interval-1]
    EqOcCflux=OcCflux[interval-1]
    EqSumOcCflux=SumOcCflux[interval-1]
    EqSedCflux=SedCflux[interval-1]
    EqSumSedCflux=SumSedCflux[interval-1]
    EqBioCflux=BioCflux[interval-1]
    SumBioCflux=sum(EqBioCflux)  
    EqCaflux_Out=Caflux_Out[interval-1]
    SumCaflux_Out=sum(EqCaflux_Out)   
    EqT1Cflux=T1Cflux[interval-1]
    EqT2Cflux=T2Cflux[interval-1]
    EqE1Cflux=E1Cflux[interval-1]
    EqE2Cflux=E2Cflux[interval-1]
    EqVolcCflux=VolcCflux[interval-1]
    EqRivCflux=RivCflux[interval-1]
    EqWeathCflux=WeathCflux[interval-1]
    EqTerrBioCflux=TerrBioCflux[interval-1]   
    EqDissC_PgC=DissCflux[interval-1]
    EqSedDissC_PgC=SedDissCflux[interval-1]
    EqCFertflux_PgC=CFertflux[interval-1]
    EqRespireflux_PgC=Respireflux[interval-1]  
    EqAnthEmit=AnthEmitOut[interval-1]
    EqAnthStock=AnthStockOut[interval-1]
    EqDeforestC=DeforestCOut[interval-1]

    # Anthropogenic emissions totals
    SumAnth_PgC=sum(AnthEmitOut)
    TotAnthInput=SumAnth_PgC
    SumDeforestC_PgC=sum(DeforestCOut)
    
    #Run rate on carbon uptake in various reservoirs
    RunTotCarbAt_PgC=TotCarbAt_PgC-TotCarbAt_PgC[0]
    RunTotCarbOc_PgC=TotCarbOc_PgC-TotCarbOc_PgC[0]
    RunCstock_PgC=Cstock_PgC-Cstock_PgC[0]
    RunSedCstock_PgC=SedCstock_PgC-SedCstock_PgC[0]
    TotCarb_PgC_delta=TotCarb_PgC[interval-1]-TotCarb_PgC[0]
    
    # Sum of other fluxes
    SumVolc=sum(VolcCflux)
    SumRivCarbCflux=sum(RivCarbCflux)
    SumCarbInputs=SumVolc+SumRivCarbCflux
    SumSedDissC=sum(SedDissCflux)
        
    # Output for "last run" files
    if StoreResults=='on':
        np.savetxt(LastDir+'LastPhos.txt',np.reshape(EqPhosOut,[7,1]))
        np.savetxt(LastDir+'LastCarb.txt',np.reshape(EqCarbOut,[7,1]))
        np.savetxt(LastDir+'LastAtCO2.txt',EqAtCO2Out)
        np.savetxt(LastDir+'LastAtO2.txt',EqAtO2Out)
        np.savetxt(LastDir+'LastAlk.txt',EqAlkOut)
        np.savetxt(LastDir+'Lastd13C.txt',np.reshape(Eqd13COut,[7,1]))
        np.savetxt(LastDir+'LastAtd13C.txt',EqAtd13COut)
        np.savetxt(LastDir+'Lastd14C.txt',np.reshape(Eqd14COut,[7,1]))
        np.savetxt(LastDir+'LastAtd14C.txt',EqAtd14COut)
        np.savetxt(LastDir+'LastFe.txt',np.reshape(EqFeOut,[7,1]))   
        np.savetxt(LastDir+'LastSi.txt',np.reshape(EqSiOut,[7,1]))
        np.savetxt(LastDir+'LastO.txt',np.reshape(EqOOut,[7,1]))
        np.savetxt(LastDir+'LastCstock1_PgC.txt',EqCstock1_PgC)
        np.savetxt(LastDir+'LastCstock2_PgC.txt',EqCstock2_PgC)
        np.savetxt(LastDir+'LastAtmosCstock_PgC.txt',EqTotCarbAt_PgC)
        np.savetxt(LastDir+'LastOceanCstock_PgC.txt',EqTotCarbOc_PgC)
        np.savetxt(LastDir+'LastCstock_PgC.txt',EqCstock_PgC)
        np.savetxt(LastDir+'LastSedCstock_PgC.txt',EqSedCstock_PgC)
        np.savetxt(LastDir+'LastContCstock_PgC.txt',EqContCstock_PgC)
        np.savetxt(LastDir+'LastRivCarbC.txt',SumRivCarbCflux)
        
    # Quick results summary
    if BatchRun=='off':    

        print('Regular outputs')
        print('Atmospheric CO2')
        print(EqAtCO2Out)
        
        print('Difference')
        print(np.round(EqAtCO2Out-(CESM[j]*1e6),2))
        print(np.round(EqAtCO2Out-(CESM[j]*1e6),3))
        
        print('Y add')
        print(y_add)

## Charts----------------------------------------------------------------------
## Remember to switch off when batch run  

# Anthropocene charts
    ls='dotted'
    mk='o'
    ms=5
    ew=1
    cap=5
    ms_data=2.0

    if AnthEmits==0:
        AnthTime=time+AnthStart # set up industrial period timing 
        if Charting=='on':
            
            

            figure1 = plt.figure(1,figsize=(8.3,11.7))
            ax1 = figure1.add_subplot(221)
            plt.title('Atmospheric carbon stock')
            ax1.plot(AnthTime,TotCarbAt_PgC)
            plt.ylabel('PgC')
            plt.xlabel('Year')
            ax2 = figure1.add_subplot(222)
            plt.title('Ocean carbon stock')
            plt.plot(AnthTime,TotCarbOc_PgC)
            plt.ylabel('PgC')
            plt.xlabel('Year')
            ax3 = figure1.add_subplot(223)
            plt.title('Terrestrial biosphere carbon stock')
            plt.plot(AnthTime,Cstock_PgC)
            plt.ylabel('PgC')
            plt.xlabel('Year')
            ax4 = figure1.add_subplot(224)
            plt.title('Ocean sediments carbon stock')
            plt.plot(AnthTime,SedCstock_PgC)
            plt.ylabel('PgC')
            plt.xlabel('Year')
            plt.tight_layout()
            plt.show()

            figure2 = plt.figure(2,figsize=(8.3,11.7))
            ax1 = figure2.add_subplot(221)
            plt.title('Net Air-sea carbon flux')
            ax1.plot(AnthTime,AtCflux)
            plt.ylabel('PgC/year')
            plt.xlabel('Year')
            ax2 = figure2.add_subplot(222)
            plt.title('Ocean sediment flux')
            plt.plot(AnthTime,SumSedCflux)
            plt.ylabel('PgC/yr')
            plt.xlabel('Year')
            ax3 = figure2.add_subplot(223)
            plt.title('Terrestrial biosphere carbon flux')
            plt.plot(AnthTime,TerrBioCflux)
            plt.plot(AnthTime,CFertflux)
            plt.plot(AnthTime,Respireflux)
            plt.ylabel('PgC/yr')
            plt.xlabel('Year')
            plt.legend(['Net','Uptake','Respiration'],loc=2)
            ax4 = figure2.add_subplot(224)
            plt.title('Other carbon fluxes')
            plt.plot(AnthTime,WeathCflux)            
            plt.plot(AnthTime,RivCflux)
            plt.plot(AnthTime,VolcCflux)
            plt.ylabel('PgC/yr')
            plt.xlabel('Year')
            plt.legend(['Weathering of silicate rocks','River flux','Volcanic'],loc=2)
            plt.tight_layout()
            plt.show()
            
            figure3 = plt.figure(3,figsize=(8.3,11.7))
            ax1 = figure3.add_subplot(321)
            plt.title('Ocean DIC')
            ax1.plot(AnthTime,CarbOut)
            ax1.errorbar(GLODAP_time,GLODAP_DIC[0],yerr=SD_GLODAP_DIC[0], ls=ls, marker=mk,markersize=ms, capsize=cap, elinewidth=ew,color='blue')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[1],yerr=SD_GLODAP_DIC[1], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='orange')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[2],yerr=SD_GLODAP_DIC[2], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='green')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[3],yerr=SD_GLODAP_DIC[3], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='red')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[4],yerr=SD_GLODAP_DIC[4], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='purple')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[5],yerr=SD_GLODAP_DIC[5], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='brown')
            ax1.errorbar(GLODAP_time,GLODAP_DIC[6],yerr=SD_GLODAP_DIC[6], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='pink')             
            plt.ylabel('umol/kg')
            plt.ylim(1500,2500)
            plt.xlabel('Year')
            ax2 = figure3.add_subplot(322)
            plt.title('Alkalinity')
            plt.plot(AnthTime,AlkOut)           
            ax2.errorbar(GLODAP_time,GLODAP_Alk[0],yerr=SD_GLODAP_Alk[0], ls=ls, marker=mk,markersize=ms, capsize=cap, elinewidth=ew,color='blue')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[1],yerr=SD_GLODAP_Alk[1], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='orange')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[2],yerr=SD_GLODAP_Alk[2], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='green')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[3],yerr=SD_GLODAP_Alk[3], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='red')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[4],yerr=SD_GLODAP_Alk[4], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='purple')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[5],yerr=SD_GLODAP_Alk[5], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='brown')
            ax2.errorbar(GLODAP_time,GLODAP_Alk[6],yerr=SD_GLODAP_Alk[6], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='pink') 
            plt.ylabel('umol/kg')
            plt.ylim(2000,2500)
            plt.xlabel('Year')
            plt.legend(['Surface eq','Northern','Intermediate','Deep','Southern','Abyssal','Subpolar'],loc=3)
            ax3 = figure3.add_subplot(323)
            plt.title('Ocean \u03B413C')
            plt.plot(AnthTime,SCarbOut)
            ax3.errorbar(GLODAP_time,GLODAP_d13C[0],yerr=SD_GLODAP_d13C[0], ls=ls, marker=mk,markersize=ms, capsize=cap, elinewidth=ew,color='blue')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[1],yerr=SD_GLODAP_d13C[1], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='orange')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[2],yerr=SD_GLODAP_d13C[2], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='green')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[3],yerr=SD_GLODAP_d13C[3], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='red')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[4],yerr=SD_GLODAP_d13C[4], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='purple')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[5],yerr=SD_GLODAP_d13C[5], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='brown')
            ax3.errorbar(GLODAP_time,GLODAP_d13C[6],yerr=SD_GLODAP_d13C[6], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='pink')       
            plt.ylabel('\u2030')
            plt.xlabel('Year')
            ax4 = figure3.add_subplot(324)
            plt.title('Ocean \u039414C')
            plt.plot(AnthTime,RCarbOut)
            ax4.errorbar(GLODAP_time,GLODAP_D14C[0],yerr=SD_GLODAP_D14C[0], ls=ls, marker=mk,markersize=ms, capsize=cap, elinewidth=ew,color='blue')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[1],yerr=SD_GLODAP_D14C[1], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='orange')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[2],yerr=SD_GLODAP_D14C[2], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='green')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[3],yerr=SD_GLODAP_D14C[3], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='red')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[4],yerr=SD_GLODAP_D14C[4], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='purple')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[5],yerr=SD_GLODAP_D14C[5], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='brown')
            ax4.errorbar(GLODAP_time,GLODAP_D14C[6],yerr=SD_GLODAP_D14C[6], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='pink') 
            plt.ylabel('\u2030')
            plt.xlabel('Year')            
            ax5 = figure3.add_subplot(325)
            plt.title('Ocean phosphate')
            ax5.plot(AnthTime,PhosOut)
            ax5.errorbar(GLODAP_time,GLODAP_Phos[0],yerr=SD_GLODAP_Phos[0], ls=ls, marker=mk,markersize=ms, capsize=cap, elinewidth=ew,color='blue')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[1],yerr=SD_GLODAP_Phos[1], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='orange')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[2],yerr=SD_GLODAP_Phos[2], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='green')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[3],yerr=SD_GLODAP_Phos[3], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='red')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[4],yerr=SD_GLODAP_Phos[4], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='purple')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[5],yerr=SD_GLODAP_Phos[5], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='brown')
            ax5.errorbar(GLODAP_time,GLODAP_Phos[6],yerr=SD_GLODAP_Phos[6], ls=ls,marker=mk,markersize=ms, capsize=cap, elinewidth=ew, color='pink') 
            plt.ylabel('umol/kg')
            plt.ylim(0.0,3.0)
            plt.xlabel('Year')
            ax6 = figure3.add_subplot(326)
            plt.title('Ocean carbonate')
            ax6.plot(AnthTime,CO23out)
            ax6.plot(GLODAP_time,GLODAP_CO23[0],'o',color='blue') 
            ax6.plot(GLODAP_time,GLODAP_CO23[1],'o',color='orange') 
            ax6.plot(GLODAP_time,GLODAP_CO23[2],'o',color='green') 
            ax6.plot(GLODAP_time,GLODAP_CO23[3],'o',color='red') 
            ax6.plot(GLODAP_time,GLODAP_CO23[4],'o',color='purple') 
            ax6.plot(GLODAP_time,GLODAP_CO23[5],'o',color='brown') 
            ax6.plot(GLODAP_time,GLODAP_CO23[6],'o',color='pink')            
            plt.ylabel('umol/kg')
            plt.ylim(0,400)
            plt.xlabel('Year')
            plt.tight_layout()            
            plt.show()
           
            figure4 = plt.figure(4,figsize=(8.3,11.7))
            ax1 = figure4.add_subplot(221)
            plt.title('Industrial emissions (incl forestry/LUC)')
            ax1.plot(AnthTime,TotalEmitOut)
            ax1.plot(AnthEmitsDat_time,AnthEmitsDat/1e3,'o',color='red',markersize=ms_data) # convert from Tg to Pg
            plt.legend(['Model','Historical Data','rcp 2.6','rcp 4.5','rcp 6.0','rcp 8.5'],loc=2)
            plt.ylabel('PgC/yr')
            plt.xlabel('Year')
            ax2 = figure4.add_subplot(222)
            plt.title('Atmosphere CO$_2$')
            plt.plot(AnthTime,AtCO2Out)
            plt.plot(AnthAtCO2_time,AnthAtCO2Dat_AtCO2,'o',color='red',markersize=ms_data)      
            plt.ylabel('ppm')
            plt.xlabel('Year')
            plt.legend(['Model','Historical Data','rcp 2.6','rcp 4.5','rcp 6.0','rcp 8.5'],loc=2)
            ax3 = figure4.add_subplot(223)
            plt.title('Atmosphere \u03B413C')
            plt.plot(AnthTime,SCAtOut)
            plt.plot(AnthAtd13C_time,AnthAtd13CDat,'o',color='red',markersize=ms_data)
            plt.ylabel('\u2030')
            plt.xlabel('Year')
            plt.legend(['Model','Historical Data'],loc=3)
            ax4 = figure4.add_subplot(224)
            plt.title('Atmosphere \u039414C')
            plt.plot(AnthTime,RCAtOut)
            plt.plot(AnthAtD14C_time,AnthAtD14CDat,'o',color='red',markersize=ms_data)
            plt.ylabel('\u2030')
            plt.xlabel('Year')
            plt.legend(['Model','Historical Data'],loc=2)
            plt.tight_layout()
            plt.show()

    if AnthEmits==1:
        if Charting=='on':
    
            figure1 = plt.figure(1,figsize=(8.3,11.7))
            ax1 = figure1.add_subplot(321)
            plt.title('Atmosphere CO$_2$', fontsize=10)
            ax1.plot(time,AtCO2Out)   
            plt.ylabel('ppm', fontsize=10)
            plt.xlabel('time (yrs)')
            plt.ylim(0,600)
            ax2 = figure1.add_subplot(323)
            plt.title('Atmosphere \u03B413C', fontsize=10)
            plt.plot(time,SCAtOut,color='green')
            plt.ylabel('\u2030', fontsize=10)
            plt.xlabel('time (years)')
            plt.ylim(-8.0,-5.0)
            ax3 = figure1.add_subplot(325)
            plt.title('Atmosphere \u039414C', fontsize=10)
            ax3.plot(time,RCAtOut,color='black')
            plt.ylabel('\u2030', fontsize=10)
            plt.xlabel('time (years)', fontsize=10)
            plt.ylim(-100,700)
            ax4 = figure1.add_subplot(322)
            plt.title('Oceanic carbonate proxy', fontsize=10)
            ax4.plot(time,CO23out)
            plt.ylabel('umol/kg', fontsize=10)
            plt.xlabel('time (years)', fontsize=10)
            ax5 = figure1.add_subplot(324)
            plt.title('Oceanic \u03B413C', fontsize=10)
            ax5.plot(time,SCarbOut)
            plt.ylabel('\u2030', fontsize=10)
            ax6 = figure1.add_subplot(326)
            plt.title('Oceanic \u039414C', fontsize=10)
            ax6.plot(time,RCarbOut)
            plt.ylabel('\u2030', fontsize=10)  
            plt.xlabel('time (years)', fontsize=10)
            plt.legend(['Surface eq','Northern','Intermediate','Deep','Southern','Abyssal','Subpolar'],loc=3)            
            plt.tight_layout()
            plt.show()    


            figure2 = plt.figure(2,figsize=(8.3,11.7))
            ax1 = figure2.add_subplot(221)
            plt.title('Ocean phosphate')
            ax1.plot(time,PhosOut)
            plt.ylabel('umol/kg')
            plt.xlabel('time (years)')
            ax2 = figure2.add_subplot(222)
            plt.title('Ocean DIC')
            plt.plot(time,CarbOut)
            plt.ylabel('umol/kg')
            plt.xlabel('time (years)')
            ax3 = figure2.add_subplot(223)
            plt.title('Ocean alkalinity')
            plt.plot(time,AlkOut)
            plt.ylabel('umol/kg')
            plt.xlabel('time (years)')
            plt.legend(['Surface eq','Northern','Intermediate','Deep','Southern','Abyssal','Subpolar'],loc=3)
            ax4 = figure2.add_subplot(224)
            plt.title('Ocean carbonate')
            plt.plot(time,CO23out)
            plt.ylabel('umol/kg')
            plt.xlabel('time (years)')
            plt.tight_layout()   
            plt.show() 

            figure3 = plt.figure(3,figsize=(8.3,11.7))
            ax1 = figure3.add_subplot(221)
            plt.title('Total carbon in Ocean and atmosphere')
            ax1.plot(time,TotCarbOc_PgCvec)
            ax1.plot(time,TotCarbAt_PgCvec)
            plt.ylabel('PgC')
            plt.xlabel('time (years)')
            plt.legend(['Ocean','Atmosphere'],loc=3)
            ax2 = figure3.add_subplot(222)
            plt.title('Total carbon in terrestrial biosphere')
            plt.plot(time,Cstock_PgC)
            plt.ylabel('PgC')
            plt.xlabel('time (years)')
            ax3 = figure3.add_subplot(223)
            plt.title('Total carbon in seds and continents')
            plt.plot(time,SedCstock_PgCvec)
            plt.plot(time,ContCstock_PgCvec)
            plt.ylabel('PgC')
            plt.xlabel('time (years)')
            plt.legend(['Sediments','Continents'],loc=3)
            ax4 = figure3.add_subplot(224)
            plt.title('Total carbon')
            plt.plot(time,TotCarb_PgCvec)
            plt.ylabel('PgC')
            plt.xlabel('time (years)')
            plt.tight_layout()
            plt.show()

            figure4 = plt.figure(4,figsize=(8.3,11.7))
            ax1 = figure4.add_subplot(221)
            plt.title('Carbon fluxes')
            ax1.plot(time, VolcCflux)  
            ax1.plot(time, WeathCflux)        
            ax1.plot(time, RivCflux)          
            ax1.plot(time, SumSedCflux)
            ax1.plot(time, VolcCflux+WeathCflux+CFertflux+Respireflux+RivCflux+SumSedCflux)        
            plt.ylim(-2,2)
            plt.ylabel('PgC/yr')
            plt.xlabel('time (years)')
            plt.legend(['Volc','Weath','River','Sed','Net total including terrestrial bios'],loc=3)          
            ax2 = figure4.add_subplot(222)
            plt.title('Terrestrial biosphere')
            plt.plot(time, CFertflux)
            plt.plot(time, Respireflux) 
            plt.ylabel('PgC/yr')
            plt.xlabel('time (years)')
            plt.legend(['Biosphere sink','Respiration'],loc=3)
            ax3 = figure4.add_subplot(224)
            plt.title('Carbonate dissolution')
            ax3.plot(time,DissCflux)
            plt.ylabel('PgC/yr')
            plt.xlabel('time (years)')
            plt.legend(['Surface eq','Northern','Intermediate','Deep','Southern','Abyssal','Subpolar'],loc=3)
            plt.tight_layout()
            plt.show()
 
    
    ## Output file-------------------------------------------------------------
    
    filename = 'output'+expname+"%04d" %nh+'.txt' 
    if BatchRun=='off':
        outfile1 = os.path.join(ResultsDir,filename)
    if BatchRun=='on':
        outfile1 = os.path.join(BatchResultsDir,filename)    
    results1=[gamma2/1e6,Psi1/1e6,Psi2/1e6,gamma1/1e6,Zed,bScale,FCAAv,PV0,PV1,PV4,PV6,kCa*secsday,n,EqAtCO2Out,EqAtd13COut,EqAtd14COut]
    results2=[EqCO23Out[0,0],EqCO23Out[0,1],EqCO23Out[0,2],EqCO23Out[0,3],EqCO23Out[0,4],EqCO23Out[0,5],EqCO23Out[0,6]]
    results3=[Eqd13COut[0,0],Eqd13COut[0,1],Eqd13COut[0,2],Eqd13COut[0,3],Eqd13COut[0,4],Eqd13COut[0,5],Eqd13COut[0,6]]     
    results4=[Eqd14COut[0,0],Eqd14COut[0,1],Eqd14COut[0,2],Eqd14COut[0,3],Eqd14COut[0,4],Eqd14COut[0,5],Eqd14COut[0,6]]
    results5=[EqPhosOut[0,0],EqPhosOut[0,1],EqPhosOut[0,2],EqPhosOut[0,3],EqPhosOut[0,4],EqPhosOut[0,5],EqPhosOut[0,6]]
    results6=[EqCarbOut[0,0],EqCarbOut[0,1],EqCarbOut[0,2],EqCarbOut[0,3],EqCarbOut[0,4],EqCarbOut[0,5],EqCarbOut[0,6]]    
    results7=[EqAlkOut[0,0],EqAlkOut[0,1],EqAlkOut[0,2],EqAlkOut[0,3],EqAlkOut[0,4],EqAlkOut[0,5],EqAlkOut[0,6]]
    results8=[EqSiOut[0,0],EqSiOut[0,1],EqSiOut[0,2],EqSiOut[0,3],EqSiOut[0,4],EqSiOut[0,5],EqSiOut[0,6]]     
    resultsT=np.concatenate((results1,results2,results3,results4,results5,results6,results7,results8),axis=0)
    results1_shape = (np.reshape(resultsT,(1,65)))
    np.savetxt(outfile1, results1_shape)

    return Parr, Carr, AtCO2, Alkarr, SCarr, SCAt, RCarr, RCAt, Fearr, Siarr, Oarr, CO23
    
    
### Allows manual model run
if BatchRun=='off':
    SCPM_run(nh=nh, expname=expname,Psi1=Psi1,Psi2=Psi2,gamma1=gamma1,gamma2=gamma2, 
             Zed=Zed,bScale=bScale,FCA0=FCA0,FCA1=FCA1,FCA4=FCA4,FCA6=FCA6,PV0=PV0,PV1=PV1,PV4=PV4,PV6=PV6,kCa=kCa,
             n=n,TerrBioCPgC=TerrBioCPgC,TC0=TC0,TC1=TC1,TC2=TC2,TC3=TC3,TC4=TC4,TC5=TC5,
             TC6=TC6, Sal0=Sal0,Sal1=Sal1,Sal4=Sal4,Sal6=Sal6,H=H,Varr=Varr,Parr=Parr,Siarr=Siarr, 
             Alkarr=Alkarr,SCarr=SCarr,RCarr=RCarr,Fearr=Fearr,Oarr=Oarr,Carr=Carr,AtCO2=AtCO2, 
             SCAt=SCAt, RCAt=RCAt,Cstock1=Cstock1,Cstock2=Cstock2,SedCstock=SedCstock,
             ContCstock=ContCstock,SedAlkstock=SedAlkstock,AnthStock=AnthStock,volcs1=volcs1,
             weathd13C=weathd13C,AtO2=AtO2)


# Ends here
