#!/usr/bin/env python3

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process

# SCPM_StartData

# This module processes data from the "GLODAP_processing.py" file and readies it
# for SCP-M model runs, typically for the baseline calibration - e.g Holocene.
# It provides corrections for the Suess effect on ocean box DIC, d13C, D14C, and also
# for bomb radiocarbon
# It is also the entry point for starting atmospheric data (references shown)

# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)


#Import modules
import math
import numpy as np

DataDir='Data_files/'


# Import GLODAP data which has been mapped into SCP-M boxes
GLODAP_Parr=np.loadtxt(DataDir+'GLODAP_Phos.txt')
GLODAP_Carr=np.loadtxt(DataDir+'GLODAP_DIC.txt')
GLODAP_Alkarr=np.loadtxt(DataDir+'GLODAP_Alk.txt')
GLODAP_SCarr=np.loadtxt(DataDir+'GLODAP_d13C.txt')
GLODAP_RCarr=np.loadtxt(DataDir+'GLODAP_d14C.txt')
GLODAP_TC=np.loadtxt(DataDir+'GLODAP_Temp.txt')
GLODAP_Sal=np.loadtxt(DataDir+'GLODAP_Sal.txt')
GLODAP_Oarr=np.loadtxt(DataDir+'GLODAP_Ox.txt')
GLODAP_Siarr=np.loadtxt(DataDir+'GLODAP_Si.txt')
GLODAP_CFC12arr=np.loadtxt(DataDir+'GLODAP_pCFC12.txt')

# Other data

# Atmospheric CO2 (ppm)
AtCO2_1990 = 352.28 # http://scrippsco2.ucsd.edu/data/atmospheric_co2/icecore_merged_products
AtCO2_1751=  276.39# 1752 http://scrippsco2.ucsd.edu/data/atmospheric_co2/icecore_merged_products
AtCO2_2016=400.66 # http://scrippsco2.ucsd.edu/data/atmospheric_co2/icecore_merged_products
#HolAvAtCO2 = 268 # Marcott et al (2014), Vostok - 11700 to 200 years ago
HolAvAtCO2=275 # Marcott et al (2014), Vostok - 6000 to 200 years ago - this is Late Holocene (being used)

#Atmospheric d13C (per mil)
Atd13C_1990 = -7.803 # Scripps/Keeling South Pole data
Atd13C_1751= -6.35 # Rubino et al (2013)
Atd13C_2016= -8.48 # Mauna Loa data
#HolAvAtd13C = -6.41 # Schmitt et al (2012) - 11700 to 200 years ago
HolAvAtd13C = -6.35 # Schmitt et al (2012) - 6000 to 200 years ago (Late Holocene)

#Atmospheric d14C (per mil)
Atd14C_1751 = 2.48 # Stuiver et al (1998)
Atd14C_1990= 151.2 # Turnbull et al(2017) 
Atd14C_1965= 641.7 # Peak value from Turnbull et al(2017) 
Atd14C_2016= 19.5 # Turnbull et al(2017)/NIWA website 
#HolAvAtd14C = 57 # Reimer Intcal 2009 - 11700 to 200 years ago
HolAvAtd14C = 20 # Reimer Intcal 2009 - 6000 to 200 years ago (Late Holocene)

#Atmospheric O2
AtO2Init=0.21

# Initial (notional) iron values in umol/kg
Fearr = np.zeros([7,1])
Fearr[0,0] = 0.5
Fearr[1,0] = 0.5
Fearr[2,0] = 0.5
Fearr[3,0] = 0.5
Fearr[4,0] = 0.5
Fearr[5,0] = 0.5
Fearr[6,0] = 0.5

# Ocean carbon corrections (Suess and weathering)
# DIC Suess Effect corrections
# Estimates from Sabine et al (2004) in umol/kg
SuessDIC=np.zeros([1,7])
SuessDIC[0,0]=45.0
SuessDIC[0,1]=45.0
SuessDIC[0,2]=25.0
SuessDIC[0,3]=2.0
SuessDIC[0,4]=5.0
SuessDIC[0,5]=0.0
SuessDIC[0,6]=45.0

HolAvDICCorr=SuessDIC
PreIndDICCorr=SuessDIC

## Suess d13C correction as per Eide et al (2017)
Eide_a=-0.001
Eide_b=-0.0004
Suessd13C=Eide_a*GLODAP_CFC12arr+Eide_b
Suessd13C=-Suessd13C # Because we're adding it back in (starts as a -ive number)
Suessd13C=np.reshape(Suessd13C,(7,1))

# D14C
SuessD14C=np.zeros([7,1])
SuessD14C[0,0]  = -275+25 # correction for bomb radiocarbon as per Key (2001), 
SuessD14C[1,0]  = -275+25 # Druffel, Toggweiler, plus Suess effect
SuessD14C[2,0]  = -275+25
SuessD14C[3,0]  = 0
SuessD14C[4,0]  = 0
SuessD14C[5,0]  = 0
SuessD14C[6,0]  = -275+25

## GLODAP start data

GLODAPStart_Parr=GLODAP_Parr
GLODAPStart_Parr=np.reshape(GLODAPStart_Parr,[7,1])
GLODAPStart_Carr=GLODAP_Carr
GLODAPStart_Carr=np.reshape(GLODAPStart_Carr,[7,1])
GLODAPStart_Alkarr=GLODAP_Alkarr
GLODAPStart_Alkarr=np.reshape(GLODAPStart_Alkarr,[7,1])
GLODAPStart_Siarr=GLODAP_Siarr
GLODAPStart_Siarr=np.reshape(GLODAPStart_Siarr,[7,1])
GLODAPStart_TC=GLODAP_TC
GLODAPStart_TC=np.reshape(GLODAPStart_TC,[7,1])
GLODAPStart_Sal=GLODAP_Sal
GLODAPStart_Sal=np.reshape(GLODAPStart_Sal,[7,1])
GLODAPStart_Oarr=GLODAP_Oarr
GLODAPStart_Oarr=np.reshape(GLODAPStart_Oarr,[7,1])
GLODAPStart_SCarr=GLODAP_SCarr
GLODAPStart_SCarr=np.reshape(GLODAP_SCarr,[7,1])
GLODAPStart_RCarr=GLODAP_RCarr
GLODAPStart_RCarr=np.reshape(GLODAP_RCarr,[7,1])
GLODAPStartSCAt = Atd13C_1990
GLODAPStartRCAt = Atd14C_1990
GLODAPStart_Fearr=Fearr

## Holocene average start data

# Corrected for Suess effect
HolAvStart_Parr=GLODAP_Parr
HolAvStart_Parr=np.reshape(HolAvStart_Parr,[7,1])
HolAvStart_Carr=GLODAP_Carr-HolAvDICCorr
HolAvStart_Carr=np.reshape(HolAvStart_Carr,[7,1])
HolAvStart_Alkarr=GLODAP_Alkarr
HolAvStart_Alkarr=np.reshape(HolAvStart_Alkarr,[7,1])
HolAvStart_Siarr=GLODAP_Siarr
HolAvStart_Siarr=np.reshape(HolAvStart_Siarr,[7,1])
HolAvStart_TC=GLODAP_TC
HolAvStart_TC=np.reshape(HolAvStart_TC,[7,1])
HolAvStart_Sal=GLODAP_Sal
HolAvStart_Sal=np.reshape(HolAvStart_Sal,[7,1])
HolAvStart_Oarr=GLODAP_Oarr
HolAvStart_Oarr=np.reshape(HolAvStart_Oarr,[7,1])
HolAvStart_SCarr = GLODAPStart_SCarr+Suessd13C
HolAvStartSCAt = HolAvAtd13C
HolAvStart_RCarr=GLODAPStart_RCarr+SuessD14C
HolAvStartRCAt = HolAvAtd14C
HolAvStart_Fearr=Fearr

## Preindustrial start data (1751)

# Corrected for Suess effect
PreIndStart_Parr=GLODAP_Parr
PreIndStart_Parr=np.reshape(PreIndStart_Parr,[7,1])
PreIndStart_Carr=GLODAP_Carr-PreIndDICCorr
PreIndStart_Carr=np.reshape(PreIndStart_Carr,[7,1])
PreIndStart_Alkarr=GLODAP_Alkarr
PreIndStart_Alkarr=np.reshape(PreIndStart_Alkarr,[7,1])
PreIndStart_Siarr=GLODAP_Siarr
PreIndStart_Siarr=np.reshape(PreIndStart_Siarr,[7,1])
PreIndStart_TC=GLODAP_TC
PreIndStart_TC=np.reshape(PreIndStart_TC,[7,1])
PreIndStart_Sal=GLODAP_Sal
PreIndStart_Sal=np.reshape(PreIndStart_Sal,[7,1])
PreIndStart_Oarr=GLODAP_Oarr
PreIndStart_Oarr=np.reshape(PreIndStart_Oarr,[7,1])
PreIndStart_SCarr = GLODAPStart_SCarr+Suessd13C
PreIndStartSCAt = Atd13C_1751
PreIndStart_RCarr=GLODAPStart_RCarr+SuessD14C
PreIndStartRCAt = Atd14C_1751
PreIndStart_Fearr=Fearr

# Datasets packaged up to use in SCP-M

GLODAPStartSet=[AtCO2_1990,Atd13C_1990,Atd14C_1990,GLODAPStart_Parr,GLODAPStart_Carr,GLODAPStart_Alkarr,GLODAPStart_SCarr,GLODAPStart_RCarr,
            GLODAPStart_Oarr,GLODAPStart_Siarr,AtO2Init]

Other_set=[Fearr]

HolAvStartSet=[HolAvAtCO2,HolAvAtd13C,HolAvAtd14C,HolAvStart_Parr,HolAvStart_Carr,HolAvStart_Alkarr,HolAvStart_SCarr,HolAvStart_RCarr,
            HolAvStart_Oarr,HolAvStart_Siarr,AtO2Init]

PreIndStartSet=[AtCO2_1751,Atd13C_1751,Atd14C_1751,PreIndStart_Parr,PreIndStart_Carr,PreIndStart_Alkarr,PreIndStart_SCarr,PreIndStart_RCarr,
            PreIndStart_Oarr,PreIndStart_Siarr,AtO2Init]




