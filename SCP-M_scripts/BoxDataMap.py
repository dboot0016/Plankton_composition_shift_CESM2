#!/usr/bin/env python3

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process

# This module maps LGM and Holocene ocean proxy data into SCP-M box space 
# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import math
import numpy as np
import csv
import pandas as pd

DataDir='Data_files/'


# Read in data
# For data sources see model documentation
Radiocarbon=pd.read_excel(DataDir+'Radiocarbon_Read.xlsx',engine="openpyxl")
StableCarbon=pd.read_excel(DataDir+'StableCarbon_Read.xlsx',engine="openpyxl")
Carbonate=pd.read_excel(DataDir+'Carbonate_Read.xlsx',engine="openpyxl")

# Age slices
LGMAge_Rad=(Radiocarbon.Age<23000) & (Radiocarbon.Age>19000)
HolAge_Rad=(Radiocarbon.Age<6000) & (Radiocarbon.Age>200)

LGMAge_Carb=(Carbonate.Age<23000) & (Carbonate.Age>19000)
HolAge_Carb=(Carbonate.Age<6000) & (Carbonate.Age>200)
# Note Peterson et al (2014 d13C data already binned into LGM and late Holocene
# slices)


#Latitude of boxes Talley (2013) interpretation (no longitude as this is zonally averaged)
Box1LatMin=-40.0
Box1LatMax=40.0
Box2LatMin=40.0
Box2LatMax=80.0 # extended for data coverage
Box3LatMin=-40.0
Box3LatMax=40.0
Box4LatMin=-61.0
Box4LatMax=60.0
Box5LatMin=-80.0
Box5LatMax=-61.0
Box6LatMin=-80.0
Box6LatMax=80.0
Box7LatMin=-61.0
Box7LatMax=-40.0

# Depth of boxes
Box1DepthMin=0.0 #Box0
Box1DepthMax=100.0 #Box0
Box2DepthMin=0.0 #Box1
Box2DepthMax=250.0 #Box1
Box3DepthMin=100.0 #Box2
Box3DepthMax=1000.0 #Box2
Box4DepthMin=1000.0 #Box3
Box4DepthMax=2500.0 #Box3
Box5DepthMin=0.0 #Box4
Box5DepthMax=2500.0 #Box4
Box6DepthMin=2500.0 #Box5
Box6DepthMax=6000.0 #Box5 # Extended depth to allow more data to be included in abyssal box (beyond average depth)
Box7DepthMin=0.0 #Box6
Box7DepthMax=250.0 #Box6

# LGM workings-----------------------------------------------------------------

# Radiocarbon

LGM_Box1Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box1LatMax) & (Radiocarbon.Lat>Box1LatMin)& (Radiocarbon.Sample_Depth>Box1DepthMin)& (Radiocarbon.Sample_Depth<Box1DepthMax)]
LGM_Box2Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box2LatMax) & (Radiocarbon.Lat>Box2LatMin)& (Radiocarbon.Sample_Depth>Box2DepthMin)& (Radiocarbon.Sample_Depth<Box2DepthMax)]
LGM_Box3Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box3LatMax) & (Radiocarbon.Lat>Box3LatMin)& (Radiocarbon.Sample_Depth>Box3DepthMin)& (Radiocarbon.Sample_Depth<Box3DepthMax)]
LGM_Box4Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box4LatMax) & (Radiocarbon.Lat>Box4LatMin)& (Radiocarbon.Sample_Depth>Box4DepthMin)& (Radiocarbon.Sample_Depth<Box4DepthMax)]
LGM_Box5Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box5LatMax) & (Radiocarbon.Lat>Box5LatMin)& (Radiocarbon.Sample_Depth>Box5DepthMin)& (Radiocarbon.Sample_Depth<Box5DepthMax)]
LGM_Box6Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box6LatMax) & (Radiocarbon.Lat>Box6LatMin)& (Radiocarbon.Sample_Depth>Box6DepthMin)& (Radiocarbon.Sample_Depth<Box6DepthMax)]
LGM_Box7Rad=Radiocarbon[LGMAge_Rad&(Radiocarbon.Lat<Box7LatMax) & (Radiocarbon.Lat>Box7LatMin)& (Radiocarbon.Sample_Depth>Box7DepthMin)& (Radiocarbon.Sample_Depth<Box7DepthMax)]

LGM_Box1Rad_Av=LGM_Box1Rad.mean()
LGM_Box2Rad_Av=LGM_Box2Rad.mean()
LGM_Box3Rad_Av=LGM_Box3Rad.mean()
LGM_Box4Rad_Av=LGM_Box4Rad.mean()
LGM_Box5Rad_Av=LGM_Box5Rad.mean()
LGM_Box6Rad_Av=LGM_Box6Rad.mean()
LGM_Box7Rad_Av=LGM_Box7Rad.mean()

LGM_Box1Rad_SD=LGM_Box1Rad.std()
LGM_Box2Rad_SD=LGM_Box2Rad.std()
LGM_Box3Rad_SD=LGM_Box3Rad.std()
LGM_Box4Rad_SD=LGM_Box4Rad.std()
LGM_Box5Rad_SD=LGM_Box5Rad.std()
LGM_Box6Rad_SD=LGM_Box6Rad.std()
LGM_Box7Rad_SD=LGM_Box7Rad.std()

# Stable carbon

LGM_Box1StableC=StableCarbon[(StableCarbon.Lat<Box1LatMax) & (StableCarbon.Lat>Box1LatMin)& (StableCarbon.Sample_Depth>Box1DepthMin)& (StableCarbon.Sample_Depth<Box1DepthMax)]
LGM_Box2StableC=StableCarbon[(StableCarbon.Lat<Box2LatMax) & (StableCarbon.Lat>Box2LatMin)& (StableCarbon.Sample_Depth>Box2DepthMin)& (StableCarbon.Sample_Depth<Box2DepthMax)]
LGM_Box3StableC=StableCarbon[(StableCarbon.Lat<Box3LatMax) & (StableCarbon.Lat>Box3LatMin)& (StableCarbon.Sample_Depth>Box3DepthMin)& (StableCarbon.Sample_Depth<Box3DepthMax)]
LGM_Box4StableC=StableCarbon[(StableCarbon.Lat<Box4LatMax) & (StableCarbon.Lat>Box4LatMin)& (StableCarbon.Sample_Depth>Box4DepthMin)& (StableCarbon.Sample_Depth<Box4DepthMax)]
LGM_Box5StableC=StableCarbon[(StableCarbon.Lat<Box5LatMax) & (StableCarbon.Lat>Box5LatMin)& (StableCarbon.Sample_Depth>Box5DepthMin)& (StableCarbon.Sample_Depth<Box5DepthMax)]
LGM_Box6StableC=StableCarbon[(StableCarbon.Lat<Box6LatMax) & (StableCarbon.Lat>Box6LatMin)& (StableCarbon.Sample_Depth>Box6DepthMin)& (StableCarbon.Sample_Depth<Box6DepthMax)]
LGM_Box7StableC=StableCarbon[(StableCarbon.Lat<Box7LatMax) & (StableCarbon.Lat>Box7LatMin)& (StableCarbon.Sample_Depth>Box7DepthMin)& (StableCarbon.Sample_Depth<Box7DepthMax)]


LGM_Box1StableC_Av=LGM_Box1StableC.mean()
LGM_Box2StableC_Av=LGM_Box2StableC.mean()
LGM_Box3StableC_Av=LGM_Box3StableC.mean()
LGM_Box4StableC_Av=LGM_Box4StableC.mean()
LGM_Box5StableC_Av=LGM_Box5StableC.mean()
LGM_Box6StableC_Av=LGM_Box6StableC.mean()
LGM_Box7StableC_Av=LGM_Box7StableC.mean()

LGM_Box1StableC_SD=LGM_Box1StableC.std()
LGM_Box2StableC_SD=LGM_Box2StableC.std()
LGM_Box3StableC_SD=LGM_Box3StableC.std()
LGM_Box4StableC_SD=LGM_Box4StableC.std()
LGM_Box5StableC_SD=LGM_Box5StableC.std()
LGM_Box6StableC_SD=LGM_Box6StableC.std()
LGM_Box7StableC_SD=LGM_Box7StableC.std()

# Carbonate

LGM_Box1Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box1LatMax) & (Carbonate.Lat>Box1LatMin)& (Carbonate.Sample_Depth>Box1DepthMin)& (Carbonate.Sample_Depth<Box1DepthMax)]
LGM_Box2Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box2LatMax) & (Carbonate.Lat>Box2LatMin)& (Carbonate.Sample_Depth>Box2DepthMin)& (Carbonate.Sample_Depth<Box2DepthMax)]
LGM_Box3Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box3LatMax) & (Carbonate.Lat>Box3LatMin)& (Carbonate.Sample_Depth>Box3DepthMin)& (Carbonate.Sample_Depth<Box3DepthMax)]
LGM_Box4Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box4LatMax) & (Carbonate.Lat>Box4LatMin)& (Carbonate.Sample_Depth>Box4DepthMin)& (Carbonate.Sample_Depth<Box4DepthMax)]
LGM_Box5Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box5LatMax) & (Carbonate.Lat>Box5LatMin)& (Carbonate.Sample_Depth>Box5DepthMin)& (Carbonate.Sample_Depth<Box5DepthMax)]
LGM_Box6Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box6LatMax) & (Carbonate.Lat>Box6LatMin)& (Carbonate.Sample_Depth>Box6DepthMin)& (Carbonate.Sample_Depth<Box6DepthMax)]
LGM_Box7Carbonate=Carbonate[LGMAge_Carb&(Carbonate.Lat<Box7LatMax) & (Carbonate.Lat>Box7LatMin)& (Carbonate.Sample_Depth>Box7DepthMin)& (Carbonate.Sample_Depth<Box7DepthMax)]

LGM_Box1Carbonate_Av=LGM_Box1Carbonate.mean()
LGM_Box2Carbonate_Av=LGM_Box2Carbonate.mean()
LGM_Box3Carbonate_Av=LGM_Box3Carbonate.mean()
LGM_Box4Carbonate_Av=LGM_Box4Carbonate.mean()
LGM_Box5Carbonate_Av=LGM_Box5Carbonate.mean()
LGM_Box6Carbonate_Av=LGM_Box6Carbonate.mean()
LGM_Box7Carbonate_Av=LGM_Box7Carbonate.mean()

LGM_Box1Carbonate_SD=LGM_Box1Carbonate.std()
LGM_Box2Carbonate_SD=LGM_Box2Carbonate.std()
LGM_Box3Carbonate_SD=LGM_Box3Carbonate.std()
LGM_Box4Carbonate_SD=LGM_Box4Carbonate.std()
LGM_Box5Carbonate_SD=LGM_Box5Carbonate.std()
LGM_Box6Carbonate_SD=LGM_Box6Carbonate.std()
LGM_Box7Carbonate_SD=LGM_Box7Carbonate.std()

# Holocene workings------------------------------------------------------------

# Radiocarbon

Hol_Box1Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box1LatMax) & (Radiocarbon.Lat>Box1LatMin)& (Radiocarbon.Sample_Depth>Box1DepthMin)& (Radiocarbon.Sample_Depth<Box1DepthMax)]
Hol_Box2Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box2LatMax) & (Radiocarbon.Lat>Box2LatMin)& (Radiocarbon.Sample_Depth>Box2DepthMin)& (Radiocarbon.Sample_Depth<Box2DepthMax)]
Hol_Box3Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box3LatMax) & (Radiocarbon.Lat>Box3LatMin)& (Radiocarbon.Sample_Depth>Box3DepthMin)& (Radiocarbon.Sample_Depth<Box3DepthMax)]
Hol_Box4Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box4LatMax) & (Radiocarbon.Lat>Box4LatMin)& (Radiocarbon.Sample_Depth>Box4DepthMin)& (Radiocarbon.Sample_Depth<Box4DepthMax)]
Hol_Box5Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box5LatMax) & (Radiocarbon.Lat>Box5LatMin)& (Radiocarbon.Sample_Depth>Box5DepthMin)& (Radiocarbon.Sample_Depth<Box5DepthMax)]
Hol_Box6Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box6LatMax) & (Radiocarbon.Lat>Box6LatMin)& (Radiocarbon.Sample_Depth>Box6DepthMin)& (Radiocarbon.Sample_Depth<Box6DepthMax)]
Hol_Box7Rad=Radiocarbon[HolAge_Rad&(Radiocarbon.Lat<Box7LatMax) & (Radiocarbon.Lat>Box7LatMin)& (Radiocarbon.Sample_Depth>Box7DepthMin)& (Radiocarbon.Sample_Depth<Box7DepthMax)]

Hol_Box1Rad_Av=Hol_Box1Rad.mean()
Hol_Box2Rad_Av=Hol_Box2Rad.mean()
Hol_Box3Rad_Av=Hol_Box3Rad.mean()
Hol_Box4Rad_Av=Hol_Box4Rad.mean()
Hol_Box5Rad_Av=Hol_Box5Rad.mean()
Hol_Box6Rad_Av=Hol_Box6Rad.mean()
Hol_Box7Rad_Av=Hol_Box7Rad.mean()

Hol_Box1Rad_SD=Hol_Box1Rad.std()
Hol_Box2Rad_SD=Hol_Box2Rad.std()
Hol_Box3Rad_SD=Hol_Box3Rad.std()
Hol_Box4Rad_SD=Hol_Box4Rad.std()
Hol_Box5Rad_SD=Hol_Box5Rad.std()
Hol_Box6Rad_SD=Hol_Box6Rad.std()
Hol_Box7Rad_SD=Hol_Box7Rad.std()

# Stable carbon

Hol_Box1StableC=StableCarbon[(StableCarbon.Lat<Box1LatMax) & (StableCarbon.Lat>Box1LatMin)& (StableCarbon.Sample_Depth>Box1DepthMin)& (StableCarbon.Sample_Depth<Box1DepthMax)]
Hol_Box2StableC=StableCarbon[(StableCarbon.Lat<Box2LatMax) & (StableCarbon.Lat>Box2LatMin)& (StableCarbon.Sample_Depth>Box2DepthMin)& (StableCarbon.Sample_Depth<Box2DepthMax)]
Hol_Box3StableC=StableCarbon[(StableCarbon.Lat<Box3LatMax) & (StableCarbon.Lat>Box3LatMin)& (StableCarbon.Sample_Depth>Box3DepthMin)& (StableCarbon.Sample_Depth<Box3DepthMax)]
Hol_Box4StableC=StableCarbon[(StableCarbon.Lat<Box4LatMax) & (StableCarbon.Lat>Box4LatMin)& (StableCarbon.Sample_Depth>Box4DepthMin)& (StableCarbon.Sample_Depth<Box4DepthMax)]
Hol_Box5StableC=StableCarbon[(StableCarbon.Lat<Box5LatMax) & (StableCarbon.Lat>Box5LatMin)& (StableCarbon.Sample_Depth>Box5DepthMin)& (StableCarbon.Sample_Depth<Box5DepthMax)]
Hol_Box6StableC=StableCarbon[(StableCarbon.Lat<Box6LatMax) & (StableCarbon.Lat>Box6LatMin)& (StableCarbon.Sample_Depth>Box6DepthMin)& (StableCarbon.Sample_Depth<Box6DepthMax)]
Hol_Box7StableC=StableCarbon[(StableCarbon.Lat<Box7LatMax) & (StableCarbon.Lat>Box7LatMin)& (StableCarbon.Sample_Depth>Box7DepthMin)& (StableCarbon.Sample_Depth<Box7DepthMax)]

Hol_Box1StableC_Av=Hol_Box1StableC.mean()
Hol_Box2StableC_Av=Hol_Box2StableC.mean()
Hol_Box3StableC_Av=Hol_Box3StableC.mean()
Hol_Box4StableC_Av=Hol_Box4StableC.mean()
Hol_Box5StableC_Av=Hol_Box5StableC.mean()
Hol_Box6StableC_Av=Hol_Box6StableC.mean()
Hol_Box7StableC_Av=Hol_Box7StableC.mean()

Hol_Box1StableC_SD=Hol_Box1StableC.std()
Hol_Box2StableC_SD=Hol_Box2StableC.std()
Hol_Box3StableC_SD=Hol_Box3StableC.std()
Hol_Box4StableC_SD=Hol_Box4StableC.std()
Hol_Box5StableC_SD=Hol_Box5StableC.std()
Hol_Box6StableC_SD=Hol_Box6StableC.std()
Hol_Box7StableC_SD=Hol_Box7StableC.std()

# Carbonate

Hol_Box1Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box1LatMax) & (Carbonate.Lat>Box1LatMin)& (Carbonate.Sample_Depth>Box1DepthMin)& (Carbonate.Sample_Depth<Box1DepthMax)]
Hol_Box2Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box2LatMax) & (Carbonate.Lat>Box2LatMin)& (Carbonate.Sample_Depth>Box2DepthMin)& (Carbonate.Sample_Depth<Box2DepthMax)]
Hol_Box3Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box3LatMax) & (Carbonate.Lat>Box3LatMin)& (Carbonate.Sample_Depth>Box3DepthMin)& (Carbonate.Sample_Depth<Box3DepthMax)]
Hol_Box4Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box4LatMax) & (Carbonate.Lat>Box4LatMin)& (Carbonate.Sample_Depth>Box4DepthMin)& (Carbonate.Sample_Depth<Box4DepthMax)]
Hol_Box5Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box5LatMax) & (Carbonate.Lat>Box5LatMin)& (Carbonate.Sample_Depth>Box5DepthMin)& (Carbonate.Sample_Depth<Box5DepthMax)]
Hol_Box6Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box6LatMax) & (Carbonate.Lat>Box6LatMin)& (Carbonate.Sample_Depth>Box6DepthMin)& (Carbonate.Sample_Depth<Box6DepthMax)]
Hol_Box7Carbonate=Carbonate[HolAge_Carb&(Carbonate.Lat<Box7LatMax) & (Carbonate.Lat>Box7LatMin)& (Carbonate.Sample_Depth>Box7DepthMin)& (Carbonate.Sample_Depth<Box7DepthMax)]

Hol_Box1Carbonate_Av=Hol_Box1Carbonate.mean()
Hol_Box2Carbonate_Av=Hol_Box2Carbonate.mean()
Hol_Box3Carbonate_Av=Hol_Box3Carbonate.mean()
Hol_Box4Carbonate_Av=Hol_Box4Carbonate.mean()
Hol_Box5Carbonate_Av=Hol_Box5Carbonate.mean()
Hol_Box6Carbonate_Av=Hol_Box6Carbonate.mean()
Hol_Box7Carbonate_Av=Hol_Box7Carbonate.mean()

Hol_Box1Carbonate_SD=Hol_Box1Carbonate.std()
Hol_Box2Carbonate_SD=Hol_Box2Carbonate.std()
Hol_Box3Carbonate_SD=Hol_Box3Carbonate.std()
Hol_Box4Carbonate_SD=Hol_Box4Carbonate.std()
Hol_Box5Carbonate_SD=Hol_Box5Carbonate.std()
Hol_Box6Carbonate_SD=Hol_Box6Carbonate.std()
Hol_Box7Carbonate_SD=Hol_Box7Carbonate.std()

# Create box arrays

LGMd14C=np.zeros([1,7])
LGMd14C[0,0]=LGM_Box1Rad_Av.D14C
LGMd14C[0,1]=LGM_Box2Rad_Av.D14C
LGMd14C[0,2]=LGM_Box3Rad_Av.D14C
LGMd14C[0,3]=LGM_Box4Rad_Av.D14C
LGMd14C[0,4]=LGM_Box5Rad_Av.D14C
LGMd14C[0,5]=LGM_Box6Rad_Av.D14C
LGMd14C[0,6]=LGM_Box7Rad_Av.D14C

Hold14C=np.zeros([1,7])
Hold14C[0,0]=Hol_Box1Rad_Av.D14C
Hold14C[0,1]=Hol_Box2Rad_Av.D14C
Hold14C[0,2]=Hol_Box3Rad_Av.D14C
Hold14C[0,3]=Hol_Box4Rad_Av.D14C
Hold14C[0,4]=Hol_Box5Rad_Av.D14C
Hold14C[0,5]=Hol_Box6Rad_Av.D14C
Hold14C[0,6]=Hol_Box7Rad_Av.D14C

LGMd13C=np.zeros([1,7])
LGMd13C[0,0]=LGM_Box1StableC_Av.LGM_d13C
LGMd13C[0,1]=LGM_Box2StableC_Av.LGM_d13C
LGMd13C[0,2]=LGM_Box3StableC_Av.LGM_d13C
LGMd13C[0,3]=LGM_Box4StableC_Av.LGM_d13C
LGMd13C[0,4]=LGM_Box5StableC_Av.LGM_d13C
LGMd13C[0,5]=LGM_Box6StableC_Av.LGM_d13C
LGMd13C[0,6]=LGM_Box7StableC_Av.LGM_d13C

Hold13C=np.zeros([1,7])
Hold13C[0,0]=Hol_Box1StableC_Av.Hol_d13C
Hold13C[0,1]=Hol_Box2StableC_Av.Hol_d13C
Hold13C[0,2]=Hol_Box3StableC_Av.Hol_d13C
Hold13C[0,3]=Hol_Box4StableC_Av.Hol_d13C
Hold13C[0,4]=Hol_Box5StableC_Av.Hol_d13C
Hold13C[0,5]=Hol_Box6StableC_Av.Hol_d13C
Hold13C[0,6]=Hol_Box7StableC_Av.Hol_d13C

LGMCO23=np.zeros([1,7])
LGMCO23[0,0]=LGM_Box1Carbonate_Av.Carbonate
LGMCO23[0,1]=LGM_Box2Carbonate_Av.Carbonate
LGMCO23[0,2]=LGM_Box3Carbonate_Av.Carbonate
LGMCO23[0,3]=LGM_Box4Carbonate_Av.Carbonate
LGMCO23[0,4]=LGM_Box5Carbonate_Av.Carbonate
LGMCO23[0,5]=LGM_Box6Carbonate_Av.Carbonate
LGMCO23[0,6]=LGM_Box7Carbonate_Av.Carbonate

HolCO23=np.zeros([1,7])
HolCO23[0,0]=Hol_Box1Carbonate_Av.Carbonate
HolCO23[0,1]=Hol_Box2Carbonate_Av.Carbonate
HolCO23[0,2]=Hol_Box3Carbonate_Av.Carbonate
HolCO23[0,3]=Hol_Box4Carbonate_Av.Carbonate
HolCO23[0,4]=Hol_Box5Carbonate_Av.Carbonate
HolCO23[0,5]=Hol_Box6Carbonate_Av.Carbonate
HolCO23[0,6]=Hol_Box7Carbonate_Av.Carbonate

# Standard deviation

SDLGMd14C=np.zeros([1,7])
SDLGMd14C[0,0]=LGM_Box1Rad_SD.D14C
SDLGMd14C[0,1]=LGM_Box2Rad_SD.D14C
SDLGMd14C[0,2]=LGM_Box3Rad_SD.D14C
SDLGMd14C[0,3]=LGM_Box4Rad_SD.D14C
SDLGMd14C[0,4]=LGM_Box5Rad_SD.D14C
SDLGMd14C[0,5]=LGM_Box6Rad_SD.D14C
SDLGMd14C[0,6]=LGM_Box7Rad_SD.D14C

SDHold14C=np.zeros([1,7])
SDHold14C[0,0]=Hol_Box1Rad_SD.D14C
SDHold14C[0,1]=Hol_Box2Rad_SD.D14C
SDHold14C[0,2]=Hol_Box3Rad_SD.D14C
SDHold14C[0,3]=Hol_Box4Rad_SD.D14C
SDHold14C[0,4]=Hol_Box5Rad_SD.D14C
SDHold14C[0,5]=Hol_Box6Rad_SD.D14C
SDHold14C[0,6]=Hol_Box7Rad_SD.D14C

SDLGMd13C=np.zeros([1,7])
SDLGMd13C[0,0]=LGM_Box1StableC_SD.LGM_d13C
SDLGMd13C[0,1]=LGM_Box2StableC_SD.LGM_d13C
SDLGMd13C[0,2]=LGM_Box3StableC_SD.LGM_d13C
SDLGMd13C[0,3]=LGM_Box4StableC_SD.LGM_d13C
SDLGMd13C[0,4]=LGM_Box5StableC_SD.LGM_d13C
SDLGMd13C[0,5]=LGM_Box6StableC_SD.LGM_d13C
SDLGMd13C[0,6]=LGM_Box7StableC_SD.LGM_d13C

SDHold13C=np.zeros([1,7])
SDHold13C[0,0]=Hol_Box1StableC_SD.Hol_d13C
SDHold13C[0,1]=Hol_Box2StableC_SD.Hol_d13C
SDHold13C[0,2]=Hol_Box3StableC_SD.Hol_d13C
SDHold13C[0,3]=Hol_Box4StableC_SD.Hol_d13C
SDHold13C[0,4]=Hol_Box5StableC_SD.Hol_d13C
SDHold13C[0,5]=Hol_Box6StableC_SD.Hol_d13C
SDHold13C[0,6]=Hol_Box7StableC_SD.Hol_d13C

SDLGMCO23=np.zeros([1,7])
SDLGMCO23[0,0]=LGM_Box1Carbonate_SD.Carbonate
SDLGMCO23[0,1]=LGM_Box2Carbonate_SD.Carbonate
SDLGMCO23[0,2]=LGM_Box3Carbonate_SD.Carbonate
SDLGMCO23[0,3]=LGM_Box4Carbonate_SD.Carbonate
SDLGMCO23[0,4]=LGM_Box5Carbonate_SD.Carbonate
SDLGMCO23[0,5]=LGM_Box6Carbonate_SD.Carbonate
SDLGMCO23[0,6]=LGM_Box7Carbonate_SD.Carbonate

SDHolCO23=np.zeros([1,7])
SDHolCO23[0,0]=Hol_Box1Carbonate_SD.Carbonate
SDHolCO23[0,1]=Hol_Box2Carbonate_SD.Carbonate
SDHolCO23[0,2]=Hol_Box3Carbonate_SD.Carbonate
SDHolCO23[0,3]=Hol_Box4Carbonate_SD.Carbonate
SDHolCO23[0,4]=Hol_Box5Carbonate_SD.Carbonate
SDHolCO23[0,5]=Hol_Box6Carbonate_SD.Carbonate
SDHolCO23[0,6]=Hol_Box7Carbonate_SD.Carbonate

## Save data files-------------------------------------------------------------

np.savetxt(DataDir+'LGMd13C_dat.txt',LGMd13C)
np.savetxt(DataDir+'Hold13C_dat.txt',Hold13C)

np.savetxt(DataDir+'LGMD14C_dat.txt',LGMd14C)
np.savetxt(DataDir+'HolD14C_dat.txt',Hold14C)

np.savetxt(DataDir+'LGMD14C_dat.txt',LGMd14C)
np.savetxt(DataDir+'HolD14C_dat.txt',Hold14C)

np.savetxt(DataDir+'SDLGMd13C_dat.txt',SDLGMd13C)
np.savetxt(DataDir+'SDHold13C_dat.txt',SDHold13C)

np.savetxt(DataDir+'SDLGMD14C_dat.txt',SDLGMd14C)
np.savetxt(DataDir+'SDHolD14C_dat.txt',SDHold14C)

## Print output----------------------------------------------------------------

#print(LGMd13C)
#print(Hold13C)
#
#print(SDLGMd13C)
#print(SDHold13C)
#
#print(LGMd14C)
#print(Hold14C)
#
#print(SDLGMd14C)
#print(SDHold14C)
#
#print(LGMCO23)
#print(HolCO23)
#
#print(SDLGMCO23)
#print(SDHolCO23)




