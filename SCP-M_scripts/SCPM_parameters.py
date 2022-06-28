#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process


## SCP-M Parameters file
# This module controls the parameter values entering the model,
# as well as the starting data fed into the model.

# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import math
import numpy as np

# Import data
from SCPM_StartData import SuessDIC
import SCPM_StartData as SCPM_Start

LastDir='Last_run/'
DataDir='Data_files/'


##-----------------CONTROLS----------------------------------------------------
# Switch to use the results from last run for starting data (must have selected option
# to store model results in the SCP-M model)
UseLastRunStart='off' # Choose 'on' or 'off'
StartData='HolAv' # Choose between 'HolAv', 'PreInd', 'GLODAP'
LGM=0 # 0 or 1 to switch on 'glacial' background settings

##-----------------------------------------------------------------------------

#Time step setup
secsyr = 60.0 * 60.0 * 24.0 * 365.0 # 31,536,000
secskyr=secsyr*1e3
secsday=60*60*24

# Box model setup
numtb=7 #number of ocean boxes

# Basic geographic setup
AREA = 3.619e14 # Surface area of the ocean https://www.ngdc.noaa.gov/mgg/global/etopo1_ocean_volumes.html
SA_Earth_m2=5.101E+14 #Surface area earth (m2)
SA_Earth_cm2=5.1E+18 #Surface area earth (cm2)
MassAt=5.1e18# NASA planet fact sheet https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
MolAt=28.97 # NASA planet fact sheet
Varrat = MassAt*1e3/MolAt # volume of the atmosphere in moles
TotVarr=1.292E+18 # Estimate of total volume of the modern ocean Levitus (1982)

#Conversions and factors
Kelv = 273.16 # convert degC into Kelvin
CVA = 1.029E-3 # convert from umol/kg to mol/m3 
CVC = 971.82 # convert from mol/m3 to umol/kg
SWD=1029 # kg/m3 seawater
M3=1/SWD
MolC=12.0107 # grams per mol of carbon
MolP=0.0322853914964# grams per mol of phosphorus
Bar=10 # Metres depth per bar of pressure
PerC=0.12 #Per cent C in CaCO3
Peta=1e15 # units of Peta
Tera=1e12 # units of Tera
Giga=1e9 # units of Giga
Avogrado=6.02E+23 # Avogrado's number atoms/mole
cm2_m2=1e4 # convert cm2 to m2

## SCP-M model parameters------------------------------------------------------
# This is where you can manually change the important flux parameters in the model
# ocean circulation (Psi1, Psi2), mixing (gamma1, gamma2) and biological pump (Zed),
# rain ratio (FCA) etc.

# Ocean physical parameters
Psi1 = 29.0*1e6 #Global overturning circuation
Psi2=19.0*1e6 # Atlantic meridional overturning circulation
gamma1 = 29.0*1e6 #Deep-abyssal vertical mixing
gamma2=40.0*1e6 #Low latitude "thermocline" mixing 

# Ocean biology parameters 
Zed=5.0 # mol C m-2 yr-1 Martin (1987) value reference at 100m depth
bScale=0.75# Martin (1987) productivity vertical decay scale
OrgDepth=100.0 #Reference depth for surface biological productiivty fluxes at 100m as per Martin (1987)

# Carbonate rain and dissolution parameters

FCA=0.07 # 0.07 #Rain ratio singular value 
COAlk = 2.0 # units of alkalinity lost per unit of CaCO3 formation in biological carbon
R = 83.144621 ##bar-cm3/(mol-K) from NIST Physical Reference Data (http://physics.nist.gov/cgi-bin/cuu/Value?r)
kCa = 0.38/secsday# 0.38 as per Sarmiento and Gruber (2005), converted from day-1 to second-1
n = 1.0 #1 as per Hales and Emerson (1997)
DissConst=2.75E-13# A CaCO3 dissolution constant parameter - tuned in model spin-up

# "Piston velocities" for air-sea gas exchange rate (m/day)
PV0=3
PV1=3
PV4=3
PV6=3

FCAadj=1.0
# Rain ratio if differentiating in each box (% of shell carbon in C org)
FCA0=0.07 #
FCA1=0.07 #
FCA4=0.07 #
FCA6=0.07 #

# Cflux scalar for each surface box, multiplicative of Cflux
# Manual exercise to set these up in the first instance, according to
# GLODAP distributions of phos, DIC and Alk
Zarr_Adj=np.zeros([7,1]) 
Zarr_Adj[0,0]=0.22
Zarr_Adj[1,0]=0.90
Zarr_Adj[4,0]=0.35
Zarr_Adj[6,0]=1.065

# Redfield ratios

rcp=130 # Organic carbon redfield ratio
rphos=1/rcp #Ratio of phosphorous to organic carbon
rcpFe = 0.1 # Iron
rcpO = -138.0# Oxygen released in photosynthesis, ratio O:P
rcN = 16.0 # Ratio of nitrogen to phosphorus as per Redfield (1963)
rcpSi = 16.0 # Assumed 1:1 with nitrate as per Brzezinski, 1985

# Some adjustment factors for undertaking manual parameter sensitivities
TC0adjust=0.0 # Add Temperature deg C low latitude surface box
TC1adjust=0.0 # Add Temperature deg C northern high latitude surface box
TC4adjust=0.0 # Add Temperature deg C Southern Ocean box
TC6adjust=0.0 # Add Temperature deg C subpolar surface box
SAadjust=1.00 # Add Multiply by ocean surface area
SalAdjust=0.00 # Add psa salinity
RCSAdjust=1.00 #Multiply by atmospheric radiocarbon production (Mariotti et al 2013 LGM best guess 1.25)
Fract=0.50 # Mulitply. Parameter to to adjust Psi1 subpolar upwelling versus southward flow into polar Southern Ocean

# Automatic application of 'LGM state' adjustments (if switched on, above) - see SCP-M documentation
if LGM == 1:
    PV4=1.00
    TC0adjust=-6.00
    TC1adjust=-6.00
    TC4adjust=-0.00
    TC6adjust=-5.50
    SAadjust=0.97
    SalAdjust=1.00
    RCSAdjust=1.25
    
# Atmospheric radiocarbon production rates
RCS_atom=1.629*RCSAdjust# atom /cm2/s #1.57 Key (2001)
RCS_bomb=700e26 # Broecker et al (1980) est bomb 14C production 1954-1963 in atoms
RCS_bombs=RCS_bomb/10/secsyr/SA_Earth_cm2 # converted to atom/cm2/s
RCS_bombs=RCS_bombs*SA_Earth_cm2/Avogrado

# CO2 Solubility parameters from Weiss (1974), in mol/kg per atm (also shown mol/l per atm)
A1_C = -60.3409#-58.0931
A2_C = 93.4517#90.5069
A3_C = 23.3585#22.2940
B1_C = 0.023517#0.027766
B2_C = -0.023656#-0.025888
B3_C = 0.0047036#0.0050578

# O2 solubility parameters from Weiss (1974), in mol/kg per atm (also shown mol/l per atm)
A1_O = -58.3877
A2_O = 85.8079
A3_O = 23.8439
B1_O = -0.034892
B2_O = 0.015568
B3_O = -0.0019387

# Carbon isotope sample reference values
Sstand=0.0112372 # PDB
Rstand=1.2e-12 #Craig (1969)

## Set up box dimensions-------------------------------------------------------

# Split of surface area of surface-exposed boxes

F0=0.75*SAadjust # Low latitude-equatorial ocean
F1=0.10*SAadjust # Northern ocean
F4=0.05*SAadjust # Southern Ocean
F5=1.00*SAadjust # Abyssal ocean volume
F6=0.10*SAadjust # Subpolar surface ocean 

# Surface area array
Sarr=np.zeros([7,1]) 
Sarr[0,0]=F0*AREA
Sarr[1,0]=F1*AREA
Sarr[2,0]=F0*AREA
Sarr[3,0]=(1-F4)*AREA
Sarr[4,0]=F4*AREA
Sarr[5,0]=F5*AREA
Sarr[6,0]=F6*AREA

SarrP=Sarr/AREA # Percentage of total surface area in each box

# Weighted average rain ratio diagnostic (using surface area splits above)
FCAAv=F0*FCA0+F1*FCA1+F4*FCA4+F6*FCA6

# Ceiling Depth array - box minimum depth.
# Cannot have zero values for the Martin scalar function
# np.ones eventually cancel out but are needed for setup.
# Used for box volume calculations and implementation of biological pump
# See model documentation for the rationale for these values (If you change them, make
# sure you do it consistently between CDepth (ceiling) and FDepth (floor))

CDepth=np.ones([7,7]) 
CDepth[0,0]=0.1e-25
CDepth[1,1]=0.1e-25
CDepth[2,0]=100
CDepth[3,0]=1000
CDepth[3,1]=250
CDepth[3,6]=250
CDepth[4,4]=0.1e-25
CDepth[5,0]=2500
CDepth[5,1]=2500
CDepth[5,4]=2500
CDepth[5,6]=2500
CDepth[6,6]=0.1e-25

#  Floor Depth array - box maximum depth
# Cannot have zero values for the Martin scalar function
# np.ones eventually cancel out but are needed for setup.
# Used for box volume calculations and implementation of biological pump
FDepth=np.ones([7,7]) 
FDepth[0,0]=100
FDepth[1,1]=250
FDepth[2,0]=1000
FDepth[3,0]=2500
FDepth[3,1]=2500
FDepth[3,6]=2500
FDepth[4,4]=2500
FDepth[5,0]=4000
FDepth[5,1]=4000
FDepth[5,4]=4000
FDepth[5,6]=4000
FDepth[6,6]=250

Thickness=(FDepth-CDepth)
MidDepth=CDepth+Thickness/2 # Median depth of each box for average pressure

# Convert box surface areas into matrix form for biological pump
# Surface area of each box with respect of overlying surface box
# determines the biological flux entering that box
SAMat=np.zeros([7,7])
SAMat[0,0]=Sarr[0,0]
SAMat[1,1]=Sarr[1,0]
SAMat[2,0]=Sarr[2,0]
SAMat[3,0]=Sarr[0,0]
SAMat[3,1]=Sarr[1,0]
SAMat[3,6]=Sarr[6,0]
SAMat[4,4]=Sarr[4,0]
SAMat[5,0]=Sarr[0,0]
SAMat[5,1]=Sarr[1,0]
SAMat[5,4]=Sarr[4,0]
SAMat[5,6]=Sarr[6,0]
SAMat[6,6]=Sarr[6,0]

SAMatP=SAMat/Sarr # Overlying surface area as percentage of total surface boxs

# Calculate mid-depth for carbonate dissolution calculations
MidDepth=MidDepth*SAMatP
MidDepth=sum(MidDepth[:,i] for i in (0,numtb-1))
MidDepth=np.reshape(MidDepth,(7,1))

#Box volumes
Varr1=SAMat*Thickness # Mulitply each box surface area by thickness
Reshape=np.ones([7,1]) # Reshape into a 7 x 1 vector
Varr=np.dot(Varr1,Reshape) # Reshape into a 7 x 1 vector
SumVarr=sum(Varr) # Total volume of the ocean

# Sediment depth and pressure for sediment carbonate dissolution
SedDepth=np.max(FDepth)
#pressure in bars at mid point depth for each box
Pbar = MidDepth/Bar 
SedPbar=FDepth[5,0]/Bar # The pressure at the bottom of the abyssal box
SedPbar=SedDepth/Bar

## Terrestrial sources and sinks-----------------------------------------------

#Initial terrestrial carbon stocks in PgC
Cstock1=96# PgC Harman et al (2011)
Cstock2=2300# PgC Harman et al (2011)

#Baseline Atmospheric CO2 for carbon fertilisation e.g. Harman et al (2011)
AtCO2baseHol=275*1e-6 # Using this
AtCO2baseLGM=195*1e-6
AtCO2baseMed=(AtCO2baseHol+AtCO2baseLGM)/2

TerrBioCPgC=40 # PgC/yr starting terrestrial biosphere flux, Harman et al (2011)
DeforestCPgC=1.5 #PgC/yr
k1A=1/3/secsyr # Starting respiration rate (fast biosphere stock), Harman et al (2011)
k2A=1/300/secsyr # Starting respiration rate (slow biosphere stock), Harman et al (2011)
B=0.8 # "beta" for terrestrial biosphere carbon fertilisation effect, Harman et al (2011)
Raupach=0.8 # Split between "fast" and "slow" terrestrial stock, # Starting respiration rate (fast biosphere stock), Harman et al (2011)

# Carbon isotope fractonation factors for terrestrial biosphere
TerrBioSC=0.975*Sstand 
TerrBioRC=0.954*Rstand

# Ocean sedimentary stocks
SedCstock=5000#PgC, surface sediments, notional value for total carbon in earth system calculations
ContCstock=5000 #PgC,stock of carbon in volcanoes, carbonate sediments. Notional value for total carbon in earth system calculations

# Volcanic emissions
volc = 6e12 # mols/annum Generally using the method of Toggweiler (2007) which is RVSIL in main model.
#For modern day scenarios, using 6e12 Marty and Tolstikhin (1998)
volcs1 = volc/secsyr; # volc in mols/sec

# Continental weathering atmospheric CO2 sink and flux into oceans
# Generally following Toggweiler (2007)
WCARB=2.0 # mol/m3/atm/yr as per Toggweiler (2007)
WCARBs=WCARB/secsyr
BSIL=0.75e-4 # mol/m3/yr as per Toggweiler (2007)
BSILs=BSIL/secsyr
WSIL=0.5 #mol/m3/atm/yr as per Toggweiler (2007)
WSILs=WSIL/secsyr

# Carbon isotops for volanic emissions and weathering
volcd13C=-6.90 #per mil
weathd13C=-6.90 #per mil
volcd13C=(volcd13C/1000+1)*Sstand
weathd13C=(weathd13C/1000+1)*Sstand
TerrRC=0.0 # 14C dead 

## Industrial emissions--------------------------------------------------------

AnthStock=0*Peta/MolC #Starting stock of industrial emissions converted from PgC into mols C

# Emissions and isotopes
AnthSC=0.978*Sstand # Fractionation factor multiplied by the standard
AnthRC=0.0*Rstand # 14C dead

## Other fluxes and parameters-------------------------------------------------

# River fluxes of phosphorus
RiverP_Tg=15.0 #15.0 Compton et al (2001) estimate at 0.7—4.8e12 g/yr reactive, 10.8-17.8e12 g/yr total
RiverP_g=RiverP_Tg*Tera
RiverP_mol=RiverP_g*MolP
RiverP_mols=RiverP_mol/secsyr

#Kwork parameters
pH = 8 # notional starting pH to allow the calculations to begin
H=10**(-pH)

#Iostope Fractionation factors

# The "thermodynamic fractionation factor" for carbon isotopes in air-sea exchange
FK = 0.9995 # 0.99915 Stable carbon as per Schmittner et al(2013)  
FKR = 0.9990 # 0.9990 Radiocarbon as per Toggweiler and Sarmiento (1985)

# Radiocarbon air-sea fractionation factors
FSAR = np.ones([7,1])
FSAR [0,0] = 0.98182
FSAR [1,0] = 0.97720
FSAR [4,0] = 0.97720
FSAR [6,0] = 0.98182

FASR = np.ones([7,1])
FASR[0,0] = 0.99786
FASR[1,0] = 0.99768
FASR[4,0] = 0.99768
FASR[6,0] = 0.99786

# Marine biological 13C/12C fractionation factor
# ref value 0.979 Schmittner (2013)
BioSCF_av = 0.9745
BioSCF = np.ones([7,1])
BioSCF[0,0] = BioSCF_av
BioSCF[1,0] = BioSCF_av
BioSCF[4,0] = BioSCF_av
BioSCF[6,0] = BioSCF_av
# The above, if it is beneficial to specify different fractionation factros

# Marine biological 14C/12C fractionation factor
# Toggweiler and Sarmiento (1985)
BioRCF=0.954

#Radiocarbon source and decay rates
RCS_mol=RCS_atom/Avogrado
RCS1At =RCS_mol*SA_Earth_cm2
RCD=(1.2097E-4/secsyr) # 1/whole life of radiocarbon which is 8267 years
RCD1 = np.ones([7,1])*RCD
RCD1At =RCD


####--------Initial element values---------------------------------------------

# Totals check on DIC Suess correction. Sabine et al (2004) 118 +/- 19
SuessDICmm=SuessDIC*CVA # convert to mol/m3
SuessDICmm=np.reshape(SuessDICmm,(7,1))
SuessMols=SuessDICmm*Varr
SuessPgC=SuessMols*MolC/Peta
SuessPgC_Tot=sum(SuessPgC)
#print(SuessPgC_Tot)

# Temperature and salinity

# Ocean box temperature in degrees Celcius
# Established outside array structure to enable batch sensitivity analysis. 
# Established in array in SCPM_model.py
# With TC adjust for scenario analysis
GLODAP_Temp=np.loadtxt(DataDir+'GLODAP_Temp.txt')
TC0=max(GLODAP_Temp[0]+TC0adjust,0) 
TC1=max(GLODAP_Temp[1]+TC1adjust,0) 
TC2=GLODAP_Temp[2]
TC3=GLODAP_Temp[3]
TC4=max(GLODAP_Temp[4]+TC4adjust,0) 
TC5=GLODAP_Temp[5]
TC6=max(GLODAP_Temp[6]+TC6adjust,0)


# Ocean box salinity
# With global SalAdjust for scenario analysis
GLODAP_Sal=np.loadtxt(DataDir+'GLODAP_Sal.txt')
Sal0=max(GLODAP_Sal[0]+SalAdjust,0)
Sal1=max(GLODAP_Sal[1]+SalAdjust,0)
Sal2=GLODAP_Sal[2]
Sal3=GLODAP_Sal[3]
Sal4=max(GLODAP_Sal[4]+SalAdjust,0)
Sal5=GLODAP_Sal[5]
Sal6=max(GLODAP_Sal[6]+SalAdjust,0)


## Data setup for SCP-M model runs---------------------------------------------

if StartData=='PreInd':
    AtCO2Init=SCPM_Start.PreIndStartSet[0]
    SCAtInit=SCPM_Start.PreIndStartSet[1]
    RCAtInit=SCPM_Start.PreIndStartSet[2]
    ParrInit=SCPM_Start.PreIndStartSet[3]
    CarrInit=SCPM_Start.PreIndStartSet[4]
    AlkarrInit=SCPM_Start.PreIndStartSet[5]
    SCarrInit=SCPM_Start.PreIndStartSet[6]
    RCarrInit=SCPM_Start.PreIndStartSet[7]
    OarrInit=SCPM_Start.PreIndStartSet[8]
    SiarrInit=SCPM_Start.PreIndStartSet[9]
    AtO2Init=SCPM_Start.PreIndStartSet[10]
    FearrInit=SCPM_Start.Other_set[0]

if StartData=='HolAv':
    AtCO2Init=SCPM_Start.HolAvStartSet[0]
    SCAtInit=SCPM_Start.HolAvStartSet[1]
    RCAtInit=SCPM_Start.HolAvStartSet[2]
    ParrInit=SCPM_Start.HolAvStartSet[3]
    CarrInit=SCPM_Start.HolAvStartSet[4]
    AlkarrInit=SCPM_Start.HolAvStartSet[5]
    SCarrInit=SCPM_Start.HolAvStartSet[6]
    RCarrInit=SCPM_Start.HolAvStartSet[7]
    OarrInit=SCPM_Start.HolAvStartSet[8]
    SiarrInit=SCPM_Start.HolAvStartSet[9]
    AtO2Init=SCPM_Start.HolAvStartSet[10]
    FearrInit=SCPM_Start.Other_set[0]

if StartData=='GLODAP':
    AtCO2Init=SCPM_Start.GLODAPStartSet[0]
    SCAtInit=SCPM_Start.GLODAPStartSet[1]
    RCAtInit=SCPM_Start.GLODAPStartSet[2]
    ParrInit=SCPM_Start.GLODAPStartSet[3]
    CarrInit=SCPM_Start.GLODAPStartSet[4]
    AlkarrInit=SCPM_Start.GLODAPStartSet[5]
    SCarrInit=SCPM_Start.GLODAPStartSet[6]
    RCarrInit=SCPM_Start.GLODAPStartSet[7]
    OarrInit=SCPM_Start.GLODAPStartSet[8]
    SiarrInit=SCPM_Start.GLODAPStartSet[9]
    AtO2Init=SCPM_Start.GLODAPStartSet[10]
    FearrInit=SCPM_Start.Other_set[0]
    AtO2Init=SCPM_Start.AtO2Init

# Setup last run data if needed
if UseLastRunStart=='on':
    LastPhos=np.loadtxt(LastDir+'LastPhos.txt')
    LastParr=np.reshape(LastPhos,[7,1])
    LastCarb=np.loadtxt(LastDir+'LastCarb.txt')
    LastCarr=np.reshape(LastCarb,[7,1])
    LastAtCO2=np.loadtxt(LastDir+'LastAtCO2.txt')
    LastAtCO2=np.reshape(LastAtCO2,[1,1])
    LastAtO2=np.loadtxt(LastDir+'LastAtO2.txt')
    LastAtO2=np.reshape(LastAtO2,[1,1])
    LastAlk=np.loadtxt(LastDir+'LastAlk.txt')
    LastAlkarr=np.reshape(LastAlk,[7,1])
    Lastd13C=np.loadtxt(LastDir+'Lastd13C.txt')
    LastSCarr=np.reshape(Lastd13C,[7,1])
    LastAtd13C=np.loadtxt(LastDir+'LastAtd13C.txt')
    LastSCAt=np.reshape(LastAtd13C,[1,1])
    Lastd14C=np.loadtxt(LastDir+'Lastd14C.txt')
    LastRCarr=np.reshape(Lastd14C,[7,1])
    LastAtd14C=np.loadtxt(LastDir+'LastAtd14C.txt')
    LastRCAt=np.reshape(LastAtd14C,[1,1])
    LastFe=np.loadtxt(LastDir+'LastFe.txt')
    LastFearr=np.reshape(LastFe,[7,1])
    LastO=np.loadtxt(LastDir+'LastO.txt')
    LastOarr=np.reshape(LastO,[7,1])    
    LastSi=np.loadtxt(LastDir+'LastSi.txt')
    LastSiarr=np.reshape(LastSi,[7,1])
    LastCstock1=np.loadtxt(LastDir+'LastCstock1_PgC.txt')
    LastCstock1=np.reshape(LastCstock1,[1,1])
    LastCstock2=np.loadtxt(LastDir+'LastCstock2_PgC.txt')
    LastCstock2=np.reshape(LastCstock2,[1,1])
    LastSedCstock=np.loadtxt(LastDir+'LastSedCstock_PgC.txt')
    LastSedCstock=np.reshape(LastSedCstock,[1,1])    
    LastContCstock=np.loadtxt(LastDir+'LastContCstock_PgC.txt')
    LastContCstock=np.reshape(LastContCstock,[1,1])       

if UseLastRunStart=='on':
    ParrInit=LastParr
    CarrInit=LastCarr
    AtCO2Init=LastAtCO2
    AtO2Init=LastAtO2
    AlkarrInit=LastAlkarr
    SCarrInit=LastSCarr
    SCAtInit=LastSCAt
    RCarrInit=LastRCarr
    RCAtInit=LastRCAt
    FearrInit=LastFearr
    SiarrInit=LastSiarr
    OarrInit=LastOarr
    Cstock1=LastCstock1
    Cstock2=LastCstock2
    SedCstock=LastSedCstock
    ContCstock=LastContCstock

# Convert input values for model runs
Parr=ParrInit*CVA # convert to mol m-3
Carr=CarrInit*CVA # convert to mol m-3
AtCO2=AtCO2Init*1e-6
AtO2=AtO2Init
Alkarr=AlkarrInit*CVA # convert to mol m-3
SCarr=(SCarrInit/1000.0+1)*Sstand*Carr
SCAt = (SCAtInit/1000+1)*Sstand*AtCO2
RCarr=(RCarrInit/1000.0+1)*Rstand*Carr
RCAt = ((RCAtInit/1000.0)+1)*Rstand*AtCO2
Fearr=FearrInit*CVA # convert to mol m-3
Siarr=SiarrInit*CVA # convert to mol m-3
Oarr=OarrInit*CVA # convert to mol m-3
Cstock1=Cstock1*Peta/MolC
Cstock2=Cstock2*Peta/MolC
SedCstock=SedCstock*Peta/MolC
ContCstock=ContCstock*Peta/MolC
SedAlkstock=SedCstock*2.0 # Not really needed