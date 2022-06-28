#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
#The [simple carbon project] model, in process

# Batch run module - SCPM_Batch

# This module is currently set up to run multi-parameter varying experiments
# It collates output into a master results file that can be analysed with data via
# the "Model.Data.py" module. 

# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)
# 

#Import modules
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import SCPM_model as SCPM
from os import listdir
from os.path import isfile, join
import warnings
warnings.simplefilter("error")


BatchResultsDir='Results/Batch/' # Results from batch run
MasterOutDir='Results/Batch_Master/' # Master batch file output


##---CONTROLS------------------------------------------------------------------

Experiment='None'
# Options ['Hol','LGM','None' (turns batch module off- leave it like this when not using,
# because you don't want to trigger a batch run when using other modules!)]

##-----------------------------------------------------------------------------

# Holocene only experiment - 4 parameters
# Really 5 parameters in the experiment but gamma2 is not varied

if Experiment=='Hol':
    
    # Holocene experiment - parameter ranges
    Psi1ex=[20.0e6,21.0e6,22.0e6,23.0e6,24.0e6,25.0e6,26.0e6,27.0e6,28.0e6,29.0e6,30.0e6,31.0e6,32.0e6,33.0e6,34.0e6,35.0e6]
    Psi2ex=[15.0e6,16.0e6,17.0e6,18.0e6,19.0e6,20.0e6,21.0e6,22.0e6,23.0e6,24.0e6,25.0e6]
    gamma1ex=[15.0e6,16.0e6,17.0e6,18.0e6,19.0e6,20.0e6,21.0e6,22.0e6,23.0e6,24.0e6,25.0e6,26.0e6,27.0e6,28.0e6,29.0e6,30.0e6]
    gamma2ex=[40.0e6]
    Zedex=[2.0,3.0,4.0,5.0,6.0,7.0]

    master='_Hol_4'
    nh = 0 # Filename counter
    for Psi1 in Psi1ex:
        for Psi2 in Psi2ex:
            for gamma1 in gamma1ex:
                for gamma2 in gamma2ex:
                    for Zed in Zedex:
                        nh += 1
                        expname='H'
                        SCPM.SCPM_run(nh,expname,Psi1, Psi2, gamma1, gamma2, 
                        Zed,bScale=SCPM.bScale,FCA0=SCPM.FCA0,FCA1=SCPM.FCA1,FCA4=SCPM.FCA4,FCA6=SCPM.FCA6,PV0=SCPM.PV0,PV1=SCPM.PV1,PV4=SCPM.PV4,PV6=SCPM.PV6,kCa=SCPM.kCa,
                        n=SCPM.n,TerrBioCPgC=SCPM.TerrBioCPgC,TC0=SCPM.TC0,TC1=SCPM.TC1,TC2=SCPM.TC2,TC3=SCPM.TC3,TC4=SCPM.TC4,TC5=SCPM.TC5,
                        TC6=SCPM.TC6,Sal0=SCPM.Sal0,Sal1=SCPM.Sal1,Sal4=SCPM.Sal4,Sal6=SCPM.Sal6,H=SCPM.H,Varr=SCPM.Varr,Parr=SCPM.Parr,Siarr=SCPM.Siarr, 
                        Alkarr=SCPM.Alkarr,SCarr=SCPM.SCarr,RCarr=SCPM.RCarr,Fearr=SCPM.Fearr,Oarr=SCPM.Oarr,Carr=SCPM.Carr,AtCO2=SCPM.AtCO2, 
                        SCAt=SCPM.SCAt,RCAt=SCPM.RCAt,Cstock1=SCPM.Cstock1,Cstock2=SCPM.Cstock2,SedCstock=SCPM.SedCstock,
                        ContCstock=SCPM.ContCstock,SedAlkstock=SCPM.SedAlkstock,AnthStock=SCPM.AnthStock,volcs1=SCPM.volcs1,
                        weathd13C=SCPM.weathd13C,AtO2=SCPM.AtO2)

    files=os.listdir(BatchResultsDir)
    with open(MasterOutDir+'Output_all'+master+'.txt', "wb") as outfile:
        for f in files:
            if f.endswith('.txt'):
                with open(os.path.join(ResultsDir, f),"rb") as infile:
                    outfile.write(infile.read())
# Once the batch experiment is completed, this script compiles all results into
# one master file, that can be be used in the "ModelData.py" optimisation.

##-----------------------------------------------------------------------------

# LGM 4 parameter
# Really 5 parameters in the experiment but gamma2 is not varied

if Experiment=='LGM':
    
    # LGM experiment - parameter ranges
    Psi1ex=[15.0e6,16.0e6,17.0e6,18.0e6,19.0e6,20.0e6,21.0e6,22.0e6,23.0e6,24.0e6,25.0e6,26.0e6,27.0e6,28.0e6,29.0e6,30.0e6]
    Psi2ex=[5.0e6,6.0e6,7.0e6,8.0e6,9.0e6,10.0e6,11.0e6,12.0e6,13.0e6,14.0e6,15.0e6,16.0e6,17.0e6,18.0e6,19.0e6,20.0e6]
    gamma1ex=[5.0e6,6.0e6,7.0e6,8.0e6,9.0e6,10.0e6,11.0e6,12.0e6,13.0e6,14.0e6,15.0e6,16.0e6,17.0e6,18.0e6,19.0e6,20.0e6,21.0e6,22.0e6,23.0e6,24.0e6,25.0e6]
    gamma2ex=[40.0e6]
    Zedex=[2.0,3.0,4.0,5.0,6.0,7.0]


    master='_LGM_4'
    nh = 0 # Filename counter
    for Psi1 in Psi1ex:
        for Psi2 in Psi2ex:
            for gamma1 in gamma1ex:
                for gamma2 in gamma2ex:
                    for Zed in Zedex:
                        nh += 1
                        expname='L'
                        SCPM.SCPM_run(nh,expname,Psi1, Psi2, gamma1, gamma2, 
                        Zed,bScale=SCPM.bScale,FCA0=SCPM.FCA0,FCA1=SCPM.FCA1,FCA4=SCPM.FCA4,FCA6=SCPM.FCA6,PV0=SCPM.PV0,PV1=SCPM.PV1,PV4=SCPM.PV4,PV6=SCPM.PV6,kCa=SCPM.kCa,
                        n=SCPM.n,TerrBioCPgC=SCPM.TerrBioCPgC,TC0=SCPM.TC0,TC1=SCPM.TC1,TC2=SCPM.TC2,TC3=SCPM.TC3,TC4=SCPM.TC4,TC5=SCPM.TC5,
                        TC6=SCPM.TC6,Sal0=SCPM.Sal0,Sal1=SCPM.Sal1,Sal4=SCPM.Sal4,Sal6=SCPM.Sal6,H=SCPM.H,Varr=SCPM.Varr,Parr=SCPM.Parr,Siarr=SCPM.Siarr, 
                        Alkarr=SCPM.Alkarr,SCarr=SCPM.SCarr,RCarr=SCPM.RCarr,Fearr=SCPM.Fearr,Oarr=SCPM.Oarr,Carr=SCPM.Carr,AtCO2=SCPM.AtCO2, 
                        SCAt=SCPM.SCAt,RCAt=SCPM.RCAt,Cstock1=SCPM.Cstock1,Cstock2=SCPM.Cstock2,SedCstock=SCPM.SedCstock,
                        ContCstock=SCPM.ContCstock,SedAlkstock=SCPM.SedAlkstock,AnthStock=SCPM.AnthStock,volcs1=SCPM.volcs1,
                        weathd13C=SCPM.weathd13C,AtO2=SCPM.AtO2)

    # Results aggregation 
    files=os.listdir(BatchResultsDir)
    with open(MasterOutDir+'Output_all'+master+'.txt', "wb") as outfile:
        for f in files:
            if f.endswith('.txt'):
                with open(os.path.join(ResultsDir, f),"rb") as infile:
                    outfile.write(infile.read())
# Once the batch experiment is completed, this script compiles all results into
# one master file, that can be be used in the "ModelData.py" optimisation.

     
