#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process

# Data management file for SCP-M
# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


from BoxDataMap import LGMd14C, Hold14C, LGMd13C, Hold13C, SDLGMd13C, SDHold13C, SDLGMd14C,SDHold14C
from BoxDataMap import LGMCO23, HolCO23, SDLGMCO23, SDHolCO23

Datadir='Data_files/'


## RCP data--------------------------------------------------------------------
# This reads in RCP data for FF, LUC emissions for input into SCP-M modern day
# simulations, as well as IPCC atmospheric CO2 projections for comparison
 
rcp_FFems=np.loadtxt(Datadir+'rcp_FFems.txt')
rcp_LUCems=np.loadtxt(Datadir+'rcp_LUCems.txt')
rcp_SST=np.loadtxt(Datadir+'rcp_SST.txt')
rcp_AtCO2=np.loadtxt(Datadir+'rcp_AtCO2.txt')
ssp_FFems=np.loadtxt(Datadir+'ssp_FFems.txt')
ssp_LUCems=np.loadtxt(Datadir+'ssp_LUCems.txt')

co2_add=np.loadtxt(Datadir+'co2_add.txt')
CESM=np.loadtxt(Datadir+'co2_cesm_esm.txt')
Psi2_CESM=np.loadtxt(Datadir+'psi2_cesm.txt')

# RCP ossil fuel emissions
rcp_time=rcp_FFems[:,0]
rcp26_FFems=rcp_FFems[:,1]
rcp45_FFems=rcp_FFems[:,2]
rcp60_FFems=rcp_FFems[:,3]
rcp85_FFems=rcp_FFems[:,4]

ssp85_FFems=ssp_FFems[:,1]
ssp85_LUCems=ssp_LUCems[:,1]

# RCP LUC emissions
rcp26_LUCems=rcp_LUCems[:,1]
rcp45_LUCems=rcp_LUCems[:,2]
rcp60_LUCems=rcp_LUCems[:,3]
rcp85_LUCems=rcp_LUCems[:,4]

# IPCC atmospheric CO2 data
rcp26_AtCO2=rcp_AtCO2[:,1]
rcp45_AtCO2=rcp_AtCO2[:,2]
rcp60_AtCO2=rcp_AtCO2[:,3]
rcp85_AtCO2=rcp_AtCO2[:,4]

# RCP sea surface temperatures
rcp26_SST=rcp_SST[:,1]
rcp45_SST=rcp_SST[:,2]
rcp60_SST=rcp_SST[:,3]
rcp85_SST=rcp_SST[:,4]

##-----------------------------------------------------------------------------
# Anthropocene historical data for atmospheric CO2, d13C and D14C
# Sources listed in SCP-M model documentation

AnthAtCO2=np.loadtxt(Datadir+'Modern_AtCO2Dat.txt')
AnthAtCO2_time=AnthAtCO2[:,0]
AnthAtCO2Dat_AtCO2=AnthAtCO2[:,1]
AnthEmitsDat=np.loadtxt(Datadir+'Hist_AtCO2Ems.txt')
AnthEmitsDat_time=AnthEmitsDat[:,0]
AnthEmitsDat=AnthEmitsDat[:,1]
AnthAtD14C=np.loadtxt(Datadir+'Modern_AtD14CDat.txt') # Includes both Turnbull (2017) NZ data and Stuiver et al (1998) data
AnthAtD14C_time=AnthAtD14C[:,0]
AnthAtD14CDat=AnthAtD14C[:,1]
AnthAtd13C=np.loadtxt(Datadir+'Modern_Atd13CDat.txt')
AnthAtd13C_time=AnthAtd13C[:,0]
AnthAtd13CDat=AnthAtd13C[:,1]


##-----------------------------------------------------------------------------
# GLODAP data - mapped into SCP-M boxes in "GLODAP_processiing.py"

GLODAP_time=1990
GLODAP_DIC=np.loadtxt(Datadir+'GLODAP_DIC.txt')
GLODAP_Alk=np.loadtxt(Datadir+'GLODAP_Alk.txt')
GLODAP_Phos=np.loadtxt(Datadir+'GLODAP_Phos.txt')
GLODAP_d13C=np.loadtxt(Datadir+'GLODAP_d13C.txt')
GLODAP_D14C=np.loadtxt(Datadir+'GLODAP_d14C.txt')
GLODAP_CO23=np.loadtxt(Datadir+'GLODAP_CO23.txt')

SD_GLODAP_DIC=np.loadtxt(Datadir+'GLODAP_DICSD.txt')
SD_GLODAP_Alk=np.loadtxt(Datadir+'GLODAP_AlkSD.txt')
SD_GLODAP_Phos=np.loadtxt(Datadir+'GLODAP_PhosSD.txt')
SD_GLODAP_d13C=np.loadtxt(Datadir+'GLODAP_d13CSD.txt')
SD_GLODAP_D14C=np.loadtxt(Datadir+'GLODAP_d14CSD.txt')

##-----------------------------------------------------------------------------
# Atmosphere data for late Holocene and LGM

# d13C - Schmitt et al (2012)
LGMAtd13C=-6.46
LGMAtd13CSD=0.01
LGMAtd13CCount=41
HolAtd13C=-6.35
HolAtd13CSD=0.01
HolAtd13CCount=56

# D14C - Reimer et al (2009)
LGMAtd14C=414
LGMAtd14CSD=32
LGMAtd14CCount=201
HolAtd14C=20.0
HolAtd14CSD=32
HolAtd14CCount=1161

# CO2 - Marcott et al (2014)
LGMAtCO2=195
LGMAtCO2SD=3
LGMAtCO2Count=31
HolAtCO2=275
HolAtCO2SD=9
HolAtCO2Count=3




