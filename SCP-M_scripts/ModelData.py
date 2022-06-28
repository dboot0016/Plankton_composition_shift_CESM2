#!/usr/bin/env python3

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process

# ModelData
# This module processes aggregated batch model results from SCP-M and then solves
# for the best fit model results and model input parameters with data in the SCP-M boxes.

# Thus module is currently enabled to evaluate Late Holocene and LGM data
# for atmospheric CO2, d13C, D14C ocean d13C, D14C and carbonate ion proxy

# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import os
import math
import numpy as np
from os import listdir
from os.path import isfile, join
import glob
import pandas as pd

InDir='Results/Batch_Master/'
OutDir='Results/ModelData_Outputs/'

## Import data 
# Atmosphere LGM and Holocene data
from DataWork import LGMAtCO2,HolAtCO2,LGMAtd13C, HolAtd13C, LGMAtd14C, HolAtd14C 
from DataWork import LGMAtCO2SD,HolAtCO2SD,LGMAtd13CSD,HolAtd13CSD,LGMAtd14CSD,HolAtd14CSD
# Ocean LGM and Holocene data
from BoxDataMap import LGMCO23, HolCO23, SDLGMCO23, SDHolCO23
from BoxDataMap import Hold14C, LGMd14C, SDLGMd14C, SDHold14C, LGMd13C,Hold13C,SDLGMd13C, SDHold13C


## CONTROLS -------------------------------------------------------------------

DataState='Hol' # Choose 'LGM' or 'Hol'
Parameter_exp='4'

# Make sure the correct "Output_all.txt" is loaded into the InDir directory
# For example, if the SCPM_batch.py module has been run correctly, for the Holocene
# 4-parameter experiment the file name should be: "Output_all_Hol_4.txt"

# Default standard deviation values for single or vary sparse 
# data observations in a box (ie no or limited use of calculated SD)
SDLGMCO23=25 # umol/kg
SDHolCO23=25 # umol/kg

##-----------------------------------------------------------------------------

# Locate model results to compare with data
AllData=np.loadtxt(InDir+'Output_all_'+DataState+'_'+Parameter_exp+'.txt')
VarAbyss=AllData[:,0]
Psi1dat=AllData[:,1]
Psi2dat=AllData[:,2]
gamma1dat=AllData[:,3]
Cfluxdat=AllData[:,4]
bScaledat=AllData[:,5]
FCAdat=AllData[:,6]
PV0dat=AllData[:,7]
PV1dat=AllData[:,8]
PV4dat=AllData[:,9]
PV6dat=AllData[:,10]
kCadat=AllData[:,11]
ndat=AllData[:,12]
AtCO2=AllData[:,13]
AtSC=AllData[:,14]
AtRC=AllData[:,15]
CO23=AllData[:,(16,17,18,19,20,21,22)]
d13C=AllData[:,(23,24,25,26,27,28,29)]
d14C=AllData[:,(30,31,32,33,34,35,36)]
Phos=AllData[:,(37,38,39,40,41,42,43)]
Carb=AllData[:,(44,45,46,47,48,49,50)]
Alk=AllData[:,(51,52,53,54,55,56,57)]
Si=AllData[:,(58,59,60,61,62,63,64)]

## Set up SCP-M model results-data residuals
LGMResidAtCO2=AtCO2-LGMAtCO2
HolResidAtCO2=AtCO2-HolAtCO2
LGMResidAtd13C=AtSC-LGMAtd13C
HolResidAtd13C=AtSC-HolAtd13C
LGMResidAtd14C=AtRC-LGMAtd14C
HolResidAtd14C=AtRC-HolAtd14C
LGMResidd13C=d13C-LGMd13C
HolResidd13C=d13C-Hold13C
LGMResidd14C=d14C-LGMd14C
HolResidd14C=d14C-Hold14C
LGMResidCO23=CO23-LGMCO23
HolResidCO23=CO23-HolCO23

## Weighting of residuals by SD
WtLGMResidAtCO2=LGMResidAtCO2/LGMAtCO2SD
WtHolResidAtCO2=HolResidAtCO2/HolAtCO2SD
WtLGMResidAtd13C=LGMResidAtd13C/LGMAtd13CSD
WtHolResidAtd13C=HolResidAtd13C/HolAtd13CSD
WtLGMResidAtd14C=LGMResidAtd14C/LGMAtd14CSD
WtHolResidAtd14C=HolResidAtd14C/HolAtd14CSD
WtLGMResidd13C=LGMResidd13C/SDLGMd13C
WtHolResidd13C=HolResidd13C/SDHold13C
WtLGMResidd14C=LGMResidd14C/SDLGMd14C
WtHolResidd14C=HolResidd14C/SDHold14C
WtLGMResidCO23=LGMResidCO23/SDLGMCO23
WtHolResidCO23=HolResidCO23/SDHolCO23

## SD-weighted squared residuals for all atmosphere variables (CO2, d13C, d14C)
TotResLGMAt=(WtLGMResidAtCO2**2+WtLGMResidAtd13C**2+WtLGMResidAtd14C**2) 
TotResHolAt=(WtHolResidAtCO2**2+WtHolResidAtd13C**2+WtHolResidAtd14C**2) 


# SD-weighted squared residuals for all ocean variables (d13C, d14C, carbonate)
# Treats nan values as zeroes to enable calculation
WtResLGMd13C=np.nansum((WtLGMResidd13C**2),axis=1)
WtResHold13C=np.nansum((WtHolResidd13C**2),axis=1)
WtResLGMd14C=np.nansum((WtLGMResidd14C**2),axis=1)
WtResHold14C=np.nansum((WtHolResidd14C**2),axis=1)
WtResLGMCO23=np.nansum((WtLGMResidCO23**2),axis=1)
WtResHolCO23=np.nansum((WtHolResidCO23**2),axis=1)

# Atmosphere and ocean boxes, sum of squared residuals, all variables
LGM=TotResLGMAt+WtResLGMd13C+WtResLGMd14C+WtResLGMCO23 
Hol=TotResHolAt+WtResHold13C+WtResHold14C+WtResHolCO23 

# Minimise sum of squared residuals for all boxes
LGM=np.where(LGM==(min(LGM)))
Hol=np.where(Hol==(min(Hol)))

# Locate the absolute minima in the data
WtLGM=AllData[LGM,:]
WtHol=AllData[Hol,:]
ResultsWtLGM = (np.reshape(WtLGM,(1,65)))
ResultsWtHol = (np.reshape(WtHol,(1,65)))

# Output of parameters and model results for the data-optimised case
if DataState=='LGM':
    np.savetxt(OutDir+'LGM_paramater_estimates.txt',ResultsWtLGM)
if DataState=='Hol':
    np.savetxt(OutDir+'Hol_parameter_estimates.txt',ResultsWtHol)

