#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process

# GLODAP_processing
# This module maps averaged GLODAP data into SCP-M box model space. 
# It also calculates box carbonate ion concentration using the GLODAP_ data.
# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import math
import numpy as np


from SCPM_parameters import Kelv, SWD, A1_C, A2_C, A3_C, B1_C, B2_C, B3_C, CDepth, FDepth
# Need these to calculate carbonate ion concentrations

DataDir='Data_files/'

# carbonate calculation inputs
pH = 8 #guess at starting pH
H=10**(-pH)

# Data downloaded from https://www.nodc.noaa.gov/ocads/oceans/GLODAPv2/
# Last accessed July 2018


# Load in and set up data
GLODAP=np.loadtxt(DataDir+'GLODAPv2.txt')
Lat=GLODAP[:,8]
Long=GLODAP[:,9]
Press=GLODAP[:,13]
Depth=GLODAP[:,14]
Temp=GLODAP[:,15]
Sal=GLODAP[:,17]
Oxygen=GLODAP[:,26]
Nitrate=GLODAP[:,31]
Nitrite=GLODAP[:,34]
Si=GLODAP[:,36]
Phos=GLODAP[:,39]
TCO2=GLODAP[:,42]
Alk=GLODAP[:,45]
CFC11=GLODAP[:,53]
pCFC12=GLODAP[:,58] 
d13C=GLODAP[:,72]
d14C=GLODAP[:,74]

# Set up duplicate finder, for intersect of depth and latitude of each box
def list_duplicates(seq):
  seen = set()
  seen_add = seen.add
  # adds all elements it doesn't know yet to seen and all other to seen_twice
  seen_twice = set( x for x in seq if x in seen or seen_add(x) )
  # turn the set into a list (as requested)
  return list(seen_twice)

# Set up SCP-M model coordinates

#Latitude SCP-M boxes (no longitude as this is zonally averaged)
Box1LatMin=-40.0
Box1LatMax=40.0
Box2LatMin=40.0
Box2LatMax=60.0 
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


# Depth SCP-M boxes
Box1DepthMin=0.0 
Box1DepthMax=100.0 
Box2DepthMin=0.0 
Box2DepthMax=250.0 
Box3DepthMin=100.0 
Box3DepthMax=1000.0 
Box4DepthMin=1000.0 
Box4DepthMax=2500.0 
Box5DepthMin=0.0 
Box5DepthMax=2500.0 
Box6DepthMin=2500.0 
Box6DepthMax=6000.0 
Box7DepthMin=0.0 
Box7DepthMax=250.0 

# Note extended depth of box 6 to allow more data to be included 
# in abyssal box (beyond average depth)

#Locate SCP-M coordinates in the data
Box1Lat=np.where(np.logical_and(Lat>=Box1LatMin, Lat<=Box1LatMax))
Box1Depth=np.where(np.logical_and(Depth>=Box1DepthMin, Depth<=Box1DepthMax))
Box2Lat=np.where(np.logical_and(Lat>=Box2LatMin, Lat<=Box2LatMax))
Box2Depth=np.where(np.logical_and(Depth>=Box2DepthMin, Depth<=Box2DepthMax))
Box3Lat=np.where(np.logical_and(Lat>=Box3LatMin, Lat<=Box3LatMax))
Box3Depth=np.where(np.logical_and(Depth>=Box3DepthMin, Depth<=Box3DepthMax))
Box4Lat=np.where(np.logical_and(Lat>=Box4LatMin, Lat<=Box4LatMax))
Box4Depth=np.where(np.logical_and(Depth>=Box4DepthMin, Depth<=Box4DepthMax))
Box5Lat=np.where(np.logical_and(Lat>=Box5LatMin, Lat<=Box5LatMax))
Box5Depth=np.where(np.logical_and(Depth>=Box5DepthMin, Depth<=Box5DepthMax))
Box6Lat=np.where(np.logical_and(Lat>=Box6LatMin, Lat<=Box6LatMax))
Box6Depth=np.where(np.logical_and(Depth>=Box6DepthMin, Depth<=Box6DepthMax))
Box7Lat=np.where(np.logical_and(Lat>=Box7LatMin, Lat<=Box7LatMax))
Box7Depth=np.where(np.logical_and(Depth>=Box7DepthMin, Depth<=Box7DepthMax))

# Find intersection of depth and latitude points
Box1=np.append(Box1Lat, Box1Depth)
Box1=list_duplicates(Box1)
Box2=np.append(Box2Lat, Box2Depth)
Box2=list_duplicates(Box2)
Box3=np.append(Box3Lat, Box3Depth)
Box3=list_duplicates(Box3)
Box4=np.append(Box4Lat, Box4Depth)
Box4=list_duplicates(Box4)
Box5=np.append(Box5Lat, Box5Depth)
Box5=list_duplicates(Box5)
Box6=np.append(Box6Lat, Box6Depth)
Box6=list_duplicates(Box6)
Box7=np.append(Box7Lat, Box7Depth)
Box7=list_duplicates(Box7)

# Elements in each box
Box1Phos=Phos[Box1]
Box2Phos=Phos[Box2]
Box3Phos=Phos[Box3]
Box4Phos=Phos[Box4]
Box5Phos=Phos[Box5]
Box6Phos=Phos[Box6]
Box7Phos=Phos[Box7]

Box1TCO2=TCO2[Box1]
Box2TCO2=TCO2[Box2]
Box3TCO2=TCO2[Box3]
Box4TCO2=TCO2[Box4]
Box5TCO2=TCO2[Box5]
Box6TCO2=TCO2[Box6]
Box7TCO2=TCO2[Box7]

Box1Alk=Alk[Box1]
Box2Alk=Alk[Box2]
Box3Alk=Alk[Box3]
Box4Alk=Alk[Box4]
Box5Alk=Alk[Box5]
Box6Alk=Alk[Box6]
Box7Alk=Alk[Box7]

Box1d14C=d14C[Box1]
Box2d14C=d14C[Box2]
Box3d14C=d14C[Box3]
Box4d14C=d14C[Box4]
Box5d14C=d14C[Box5]
Box6d14C=d14C[Box6]
Box7d14C=d14C[Box7]

Box1d13C=d13C[Box1]
Box2d13C=d13C[Box2]
Box3d13C=d13C[Box3]
Box4d13C=d13C[Box4]
Box5d13C=d13C[Box5]
Box6d13C=d13C[Box6]
Box7d13C=d13C[Box7]

Box1Temp=Temp[Box1]
Box2Temp=Temp[Box2]
Box3Temp=Temp[Box3]
Box4Temp=Temp[Box4]
Box5Temp=Temp[Box5]
Box6Temp=Temp[Box6]
Box7Temp=Temp[Box7]

Box1Sal=Sal[Box1]
Box2Sal=Sal[Box2]
Box3Sal=Sal[Box3]
Box4Sal=Sal[Box4]
Box5Sal=Sal[Box5]
Box6Sal=Sal[Box6]
Box7Sal=Sal[Box7]

Box1Press=Press[Box1]
Box2Press=Press[Box2]
Box3Press=Press[Box3]
Box4Press=Press[Box4]
Box5Press=Press[Box5]
Box6Press=Press[Box6]
Box7Press=Press[Box7]

Box1Ox=Oxygen[Box1]
Box2Ox=Oxygen[Box2]
Box3Ox=Oxygen[Box3]
Box4Ox=Oxygen[Box4]
Box5Ox=Oxygen[Box5]
Box6Ox=Oxygen[Box6]
Box7Ox=Oxygen[Box7]

Box1Si=Si[Box1]
Box2Si=Si[Box2]
Box3Si=Si[Box3]
Box4Si=Si[Box4]
Box5Si=Si[Box5]
Box6Si=Si[Box6]
Box7Si=Si[Box7]

Box1pCFC12=pCFC12[Box1]
Box2pCFC12=pCFC12[Box2]
Box3pCFC12=pCFC12[Box3]
Box4pCFC12=pCFC12[Box4]
Box5pCFC12=pCFC12[Box5]
Box6pCFC12=pCFC12[Box6]
Box7pCFC12=pCFC12[Box7]

# Subset of values (net of nulls or negatives - they turn up as negative values, 
# or -9 and -999 in the case of the isotopes)

Box1PhosPos=[x for x in Box1Phos if x >= 0.0 ]
Box2PhosPos=[x for x in Box2Phos if x >= 0.0 ]
Box3PhosPos=[x for x in Box3Phos if x >= 0.0 ]
Box4PhosPos=[x for x in Box4Phos if x >= 0.0 ]
Box5PhosPos=[x for x in Box5Phos if x >= 0.0 ]
Box6PhosPos=[x for x in Box6Phos if x >= 0.0 ]
Box7PhosPos=[x for x in Box7Phos if x >= 0.0 ]

Box1TCO2Pos=[x for x in Box1TCO2 if x >= 0.0 ]
Box2TCO2Pos=[x for x in Box2TCO2 if x >= 0.0 ]
Box3TCO2Pos=[x for x in Box3TCO2 if x >= 0.0 ]
Box4TCO2Pos=[x for x in Box4TCO2 if x >= 0.0 ]
Box5TCO2Pos=[x for x in Box5TCO2 if x >= 0.0 ]
Box6TCO2Pos=[x for x in Box6TCO2 if x >= 0.0 ]
Box7TCO2Pos=[x for x in Box7TCO2 if x >= 0.0 ]

Box1AlkPos=[x for x in Box1Alk if x >= 0.0 ]
Box2AlkPos=[x for x in Box2Alk if x >= 0.0 ]
Box3AlkPos=[x for x in Box3Alk if x >= 0.0 ]
Box4AlkPos=[x for x in Box4Alk if x >= 0.0 ]
Box5AlkPos=[x for x in Box5Alk if x >= 0.0 ]
Box6AlkPos=[x for x in Box6Alk if x >= 0.0 ]
Box7AlkPos=[x for x in Box7Alk if x >= 0.0 ]

Box1d14CPos=[x for x in Box1d14C if x > -999.0]
Box2d14CPos=[x for x in Box2d14C if x > -999.0]
Box3d14CPos=[x for x in Box3d14C if x > -999.0]
Box4d14CPos=[x for x in Box4d14C if x > -999.0]
Box5d14CPos=[x for x in Box5d14C if x > -999.0]
Box6d14CPos=[x for x in Box6d14C if x > -999.0]
Box7d14CPos=[x for x in Box7d14C if x > -999.0]

Box1d13CPos=[x for x in Box1d13C if x >-9.0]
Box2d13CPos=[x for x in Box2d13C if x >-9.0]
Box3d13CPos=[x for x in Box3d13C if x >-9.0]
Box4d13CPos=[x for x in Box4d13C if x >-9.0]
Box5d13CPos=[x for x in Box5d13C if x >-9.0]
Box6d13CPos=[x for x in Box6d13C if x >-9.0]
Box7d13CPos=[x for x in Box7d13C if x >-9.0]

Box1TempPos=[x for x in Box1Temp if x >0.0]
Box2TempPos=[x for x in Box2Temp if x >0.0]
Box3TempPos=[x for x in Box3Temp if x >0.0]
Box4TempPos=[x for x in Box4Temp if x >0.0]
Box5TempPos=[x for x in Box5Temp if x >0.0]
Box6TempPos=[x for x in Box6Temp if x >0.0]
Box7TempPos=[x for x in Box7Temp if x >0.0]

Box1SalPos=[x for x in Box1Sal if x >0.0]
Box2SalPos=[x for x in Box2Sal if x >0.0]
Box3SalPos=[x for x in Box3Sal if x >0.0]
Box4SalPos=[x for x in Box4Sal if x >0.0]
Box5SalPos=[x for x in Box5Sal if x >0.0]
Box6SalPos=[x for x in Box6Sal if x >0.0]
Box7SalPos=[x for x in Box7Sal if x >0.0]

Box1PressPos=[x for x in Box1Press if x >0.0]
Box2PressPos=[x for x in Box2Press if x >0.0]
Box3PressPos=[x for x in Box3Press if x >0.0]
Box4PressPos=[x for x in Box4Press if x >0.0]
Box5PressPos=[x for x in Box5Press if x >0.0]
Box6PressPos=[x for x in Box6Press if x >0.0]
Box7PressPos=[x for x in Box7Press if x >0.0]

Box1OxPos=[x for x in Box1Ox if x >0.0]
Box2OxPos=[x for x in Box2Ox if x >0.0]
Box3OxPos=[x for x in Box3Ox if x >0.0]
Box4OxPos=[x for x in Box4Ox if x >0.0]
Box5OxPos=[x for x in Box5Ox if x >0.0]
Box6OxPos=[x for x in Box6Ox if x >0.0]
Box7OxPos=[x for x in Box7Ox if x >0.0]

Box1SiPos=[x for x in Box1Si if x >0.0]
Box2SiPos=[x for x in Box2Si if x >0.0]
Box3SiPos=[x for x in Box3Si if x >0.0]
Box4SiPos=[x for x in Box4Si if x >0.0]
Box5SiPos=[x for x in Box5Si if x >0.0]
Box6SiPos=[x for x in Box6Si if x >0.0]
Box7SiPos=[x for x in Box7Si if x >0.0]

Box1pCFC12Pos=[x for x in Box1pCFC12 if x >0.0]
Box2pCFC12Pos=[x for x in Box2pCFC12 if x >0.0]
Box3pCFC12Pos=[x for x in Box3pCFC12 if x >0.0]
Box4pCFC12Pos=[x for x in Box4pCFC12 if x >0.0]
Box5pCFC12Pos=[x for x in Box5pCFC12 if x >0.0]
Box6pCFC12Pos=[x for x in Box6pCFC12 if x >0.0]
Box7pCFC12Pos=[x for x in Box7pCFC12 if x >0.0]

# Simple zonal average of element concentration in each box--------------------

Box1PhosAv=np.average(Box1PhosPos)
Box2PhosAv=np.average(Box2PhosPos)
Box3PhosAv=np.average(Box3PhosPos)
Box4PhosAv=np.average(Box4PhosPos)
Box5PhosAv=np.average(Box5PhosPos)
Box6PhosAv=np.average(Box6PhosPos)
Box7PhosAv=np.average(Box7PhosPos)

Box1TCO2Av=np.average(Box1TCO2Pos)
Box2TCO2Av=np.average(Box2TCO2Pos)
Box3TCO2Av=np.average(Box3TCO2Pos)
Box4TCO2Av=np.average(Box4TCO2Pos)
Box5TCO2Av=np.average(Box5TCO2Pos)
Box6TCO2Av=np.average(Box6TCO2Pos)
Box7TCO2Av=np.average(Box7TCO2Pos)

Box1AlkAv=np.average(Box1AlkPos)
Box2AlkAv=np.average(Box2AlkPos)
Box3AlkAv=np.average(Box3AlkPos)
Box4AlkAv=np.average(Box4AlkPos)
Box5AlkAv=np.average(Box5AlkPos)
Box6AlkAv=np.average(Box6AlkPos)
Box7AlkAv=np.average(Box7AlkPos)

Box1d14CAv=np.average(Box1d14CPos)
Box2d14CAv=np.average(Box2d14CPos)
Box3d14CAv=np.average(Box3d14CPos)
Box4d14CAv=np.average(Box4d14CPos)
Box5d14CAv=np.average(Box5d14CPos)
Box6d14CAv=np.average(Box6d14CPos)
Box7d14CAv=np.average(Box7d14CPos)

Box1d13CAv=np.average(Box1d13CPos)
Box2d13CAv=np.average(Box2d13CPos)
Box3d13CAv=np.average(Box3d13CPos)
Box4d13CAv=np.average(Box4d13CPos)
Box5d13CAv=np.average(Box5d13CPos)
Box6d13CAv=np.average(Box6d13CPos)
Box7d13CAv=np.average(Box7d13CPos)

Box1TempAv=np.average(Box1TempPos)
Box2TempAv=np.average(Box2TempPos)
Box3TempAv=np.average(Box3TempPos)
Box4TempAv=np.average(Box4TempPos)
Box5TempAv=np.average(Box5TempPos)
Box6TempAv=np.average(Box6TempPos)
Box7TempAv=np.average(Box7TempPos)

Box1SalAv=np.average(Box1SalPos)
Box2SalAv=np.average(Box2SalPos)
Box3SalAv=np.average(Box3SalPos)
Box4SalAv=np.average(Box4SalPos)
Box5SalAv=np.average(Box5SalPos)
Box6SalAv=np.average(Box6SalPos)
Box7SalAv=np.average(Box7SalPos)

Box1PressAv=np.average(Box1PressPos)
Box2PressAv=np.average(Box2PressPos)
Box3PressAv=np.average(Box3PressPos)
Box4PressAv=np.average(Box4PressPos)
Box5PressAv=np.average(Box5PressPos)
Box6PressAv=np.average(Box6PressPos)
Box7PressAv=np.average(Box7PressPos)

Box1OxAv=np.average(Box1OxPos)
Box2OxAv=np.average(Box2OxPos)
Box3OxAv=np.average(Box3OxPos)
Box4OxAv=np.average(Box4OxPos)
Box5OxAv=np.average(Box5OxPos)
Box6OxAv=np.average(Box6OxPos)
Box7OxAv=np.average(Box7OxPos)

Box1SiAv=np.average(Box1SiPos)
Box2SiAv=np.average(Box2SiPos)
Box3SiAv=np.average(Box3SiPos)
Box4SiAv=np.average(Box4SiPos)
Box5SiAv=np.average(Box5SiPos)
Box6SiAv=np.average(Box6SiPos)
Box7SiAv=np.average(Box7SiPos)

Box1pCFC12Av=np.average(Box1pCFC12Pos)
Box2pCFC12Av=np.average(Box2pCFC12Pos)
Box3pCFC12Av=np.average(Box3pCFC12Pos)
Box4pCFC12Av=np.average(Box4pCFC12Pos)
Box5pCFC12Av=np.average(Box5pCFC12Pos)
Box6pCFC12Av=np.average(Box6pCFC12Pos)
Box7pCFC12Av=np.average(Box7pCFC12Pos)

# Standard deviations

Box1PhosSD=np.std(Box1PhosPos)
Box2PhosSD=np.std(Box2PhosPos)
Box3PhosSD=np.std(Box3PhosPos)
Box4PhosSD=np.std(Box4PhosPos)
Box5PhosSD=np.std(Box5PhosPos)
Box6PhosSD=np.std(Box6PhosPos)
Box7PhosSD=np.std(Box7PhosPos)

Box1TCO2SD=np.std(Box1TCO2Pos)
Box2TCO2SD=np.std(Box2TCO2Pos)
Box3TCO2SD=np.std(Box3TCO2Pos)
Box4TCO2SD=np.std(Box4TCO2Pos)
Box5TCO2SD=np.std(Box5TCO2Pos)
Box6TCO2SD=np.std(Box6TCO2Pos)
Box7TCO2SD=np.std(Box7TCO2Pos)

Box1AlkSD=np.std(Box1AlkPos)
Box2AlkSD=np.std(Box2AlkPos)
Box3AlkSD=np.std(Box3AlkPos)
Box4AlkSD=np.std(Box4AlkPos)
Box5AlkSD=np.std(Box5AlkPos)
Box6AlkSD=np.std(Box6AlkPos)
Box7AlkSD=np.std(Box7AlkPos)

Box1d14CSD=np.std(Box1d14CPos)
Box2d14CSD=np.std(Box2d14CPos)
Box3d14CSD=np.std(Box3d14CPos)
Box4d14CSD=np.std(Box4d14CPos)
Box5d14CSD=np.std(Box5d14CPos)
Box6d14CSD=np.std(Box6d14CPos)
Box7d14CSD=np.std(Box7d14CPos)

Box1d13CSD=np.std(Box1d13CPos)
Box2d13CSD=np.std(Box2d13CPos)
Box3d13CSD=np.std(Box3d13CPos)
Box4d13CSD=np.std(Box4d13CPos)
Box5d13CSD=np.std(Box5d13CPos)
Box6d13CSD=np.std(Box6d13CPos)
Box7d13CSD=np.std(Box7d13CPos)

Box1TempSD=np.std(Box1TempPos)
Box2TempSD=np.std(Box2TempPos)
Box3TempSD=np.std(Box3TempPos)
Box4TempSD=np.std(Box4TempPos)
Box5TempSD=np.std(Box5TempPos)
Box6TempSD=np.std(Box6TempPos)
Box7TempSD=np.std(Box7TempPos)

Box1SalSD=np.std(Box1SalPos)
Box2SalSD=np.std(Box2SalPos)
Box3SalSD=np.std(Box3SalPos)
Box4SalSD=np.std(Box4SalPos)
Box5SalSD=np.std(Box5SalPos)
Box6SalSD=np.std(Box6SalPos)
Box7SalSD=np.std(Box7SalPos)

Box1PressSD=np.std(Box1PressPos)
Box2PressSD=np.std(Box2PressPos)
Box3PressSD=np.std(Box3PressPos)
Box4PressSD=np.std(Box4PressPos)
Box5PressSD=np.std(Box5PressPos)
Box6PressSD=np.std(Box6PressPos)
Box7PressSD=np.std(Box7PressPos)

Box1OxSD=np.std(Box1OxPos)
Box2OxSD=np.std(Box2OxPos)
Box3OxSD=np.std(Box3OxPos)
Box4OxSD=np.std(Box4OxPos)
Box5OxSD=np.std(Box5OxPos)
Box6OxSD=np.std(Box6OxPos)
Box7OxSD=np.std(Box7OxPos)

Box1SiSD=np.std(Box1SiPos)
Box2SiSD=np.std(Box2SiPos)
Box3SiSD=np.std(Box3SiPos)
Box4SiSD=np.std(Box4SiPos)
Box5SiSD=np.std(Box5SiPos)
Box6SiSD=np.std(Box6SiPos)
Box7SiSD=np.std(Box7SiPos)

Box1pCFC12SD=np.std(Box1pCFC12Pos)
Box2pCFC12SD=np.std(Box2pCFC12Pos)
Box3pCFC12SD=np.std(Box3pCFC12Pos)
Box4pCFC12SD=np.std(Box4pCFC12Pos)
Box5pCFC12SD=np.std(Box5pCFC12Pos)
Box6pCFC12SD=np.std(Box6pCFC12Pos)
Box7pCFC12SD=np.std(Box7pCFC12Pos)

#Set up arrays-----------------------------------------------------------------

## Box averages
Phos=np.zeros([1,7])
Phos[0,0]=Box1PhosAv
Phos[0,1]=Box2PhosAv
Phos[0,2]=Box3PhosAv
Phos[0,3]=Box4PhosAv
Phos[0,4]=Box5PhosAv
Phos[0,5]=Box6PhosAv
Phos[0,6]=Box7PhosAv
Phos=np.round(Phos,2)

TCO2=np.zeros([1,7])
TCO2[0,0]=Box1TCO2Av
TCO2[0,1]=Box2TCO2Av
TCO2[0,2]=Box3TCO2Av
TCO2[0,3]=Box4TCO2Av
TCO2[0,4]=Box5TCO2Av
TCO2[0,5]=Box6TCO2Av
TCO2[0,6]=Box7TCO2Av
TCO2=np.round(TCO2,2)

Alk=np.zeros([1,7])
Alk[0,0]=Box1AlkAv
Alk[0,1]=Box2AlkAv
Alk[0,2]=Box3AlkAv
Alk[0,3]=Box4AlkAv
Alk[0,4]=Box5AlkAv
Alk[0,5]=Box6AlkAv
Alk[0,6]=Box7AlkAv
Alk=np.round(Alk,2)

d14C=np.zeros([1,7])
d14C[0,0]=Box1d14CAv
d14C[0,1]=Box2d14CAv
d14C[0,2]=Box3d14CAv
d14C[0,3]=Box4d14CAv
d14C[0,4]=Box5d14CAv
d14C[0,5]=Box6d14CAv
d14C[0,6]=Box7d14CAv
d14C=np.round(d14C,2)

d13C=np.zeros([1,7])
d13C[0,0]=Box1d13CAv
d13C[0,1]=Box2d13CAv
d13C[0,2]=Box3d13CAv
d13C[0,3]=Box4d13CAv
d13C[0,4]=Box5d13CAv
d13C[0,5]=Box6d13CAv
d13C[0,6]=Box7d13CAv
d13C=np.round(d13C,2)

Temp=np.zeros([1,7])
Temp[0,0]=Box1TempAv
Temp[0,1]=Box2TempAv
Temp[0,2]=Box3TempAv
Temp[0,3]=Box4TempAv
Temp[0,4]=Box5TempAv
Temp[0,5]=Box6TempAv
Temp[0,6]=Box7TempAv
Temp=np.round(Temp,2)

Sal=np.zeros([1,7])
Sal[0,0]=Box1SalAv
Sal[0,1]=Box2SalAv
Sal[0,2]=Box3SalAv
Sal[0,3]=Box4SalAv
Sal[0,4]=Box5SalAv
Sal[0,5]=Box6SalAv
Sal[0,6]=Box7SalAv
Sal=np.round(Sal,2)

Press=np.zeros([1,7])
Press[0,0]=Box1PressAv
Press[0,1]=Box2PressAv
Press[0,2]=Box3PressAv
Press[0,3]=Box4PressAv
Press[0,4]=Box5PressAv
Press[0,5]=Box6PressAv
Press[0,6]=Box7PressAv
Press=np.round(Press,2)

Ox=np.zeros([1,7])
Ox[0,0]=Box1OxAv
Ox[0,1]=Box2OxAv
Ox[0,2]=Box3OxAv
Ox[0,3]=Box4OxAv
Ox[0,4]=Box5OxAv
Ox[0,5]=Box6OxAv
Ox[0,6]=Box7OxAv
Ox=np.round(Ox,1)

Si=np.zeros([1,7])
Si[0,0]=Box1SiAv
Si[0,1]=Box2SiAv
Si[0,2]=Box3SiAv
Si[0,3]=Box4SiAv
Si[0,4]=Box5SiAv
Si[0,5]=Box6SiAv
Si[0,6]=Box7SiAv
Si=np.round(Si,1)

pCFC12=np.zeros([1,7])
pCFC12[0,0]=Box1pCFC12Av
pCFC12[0,1]=Box2pCFC12Av
pCFC12[0,2]=Box3pCFC12Av
pCFC12[0,3]=Box4pCFC12Av
pCFC12[0,4]=Box5pCFC12Av
pCFC12[0,5]=Box6pCFC12Av
pCFC12[0,6]=Box7pCFC12Av
pCFC12=np.round(pCFC12,6)

# Standard deviations
SDPhos=np.zeros([1,7])
SDPhos[0,0]=Box1PhosSD
SDPhos[0,1]=Box2PhosSD
SDPhos[0,2]=Box3PhosSD
SDPhos[0,3]=Box4PhosSD
SDPhos[0,4]=Box5PhosSD
SDPhos[0,5]=Box6PhosSD
SDPhos[0,6]=Box7PhosSD
SDPhos=np.round(SDPhos,2)

SDTCO2=np.zeros([1,7])
SDTCO2[0,0]=Box1TCO2SD
SDTCO2[0,1]=Box2TCO2SD
SDTCO2[0,2]=Box3TCO2SD
SDTCO2[0,3]=Box4TCO2SD
SDTCO2[0,4]=Box5TCO2SD
SDTCO2[0,5]=Box6TCO2SD
SDTCO2[0,6]=Box7TCO2SD
SDTCO2=np.round(SDTCO2,2)

SDAlk=np.zeros([1,7])
SDAlk[0,0]=Box1AlkSD
SDAlk[0,1]=Box2AlkSD
SDAlk[0,2]=Box3AlkSD
SDAlk[0,3]=Box4AlkSD
SDAlk[0,4]=Box5AlkSD
SDAlk[0,5]=Box6AlkSD
SDAlk[0,6]=Box7AlkSD
SDAlk=np.round(SDAlk,2)

SDd14C=np.zeros([1,7])
SDd14C[0,0]=Box1d14CSD
SDd14C[0,1]=Box2d14CSD
SDd14C[0,2]=Box3d14CSD
SDd14C[0,3]=Box4d14CSD
SDd14C[0,4]=Box5d14CSD
SDd14C[0,5]=Box6d14CSD
SDd14C[0,6]=Box7d14CSD
SDd14C=np.round(SDd14C,2)

SDd13C=np.zeros([1,7])
SDd13C[0,0]=Box1d13CSD
SDd13C[0,1]=Box2d13CSD
SDd13C[0,2]=Box3d13CSD
SDd13C[0,3]=Box4d13CSD
SDd13C[0,4]=Box5d13CSD
SDd13C[0,5]=Box6d13CSD
SDd13C[0,6]=Box7d13CSD
SDd13C=np.round(SDd13C,2)

SDTemp=np.zeros([1,7])
SDTemp[0,0]=Box1TempSD
SDTemp[0,1]=Box2TempSD
SDTemp[0,2]=Box3TempSD
SDTemp[0,3]=Box4TempSD
SDTemp[0,4]=Box5TempSD
SDTemp[0,5]=Box6TempSD
SDTemp[0,6]=Box7TempSD
SDTemp=np.round(SDTemp,2)

SDSal=np.zeros([1,7])
SDSal[0,0]=Box1SalSD
SDSal[0,1]=Box2SalSD
SDSal[0,2]=Box3SalSD
SDSal[0,3]=Box4SalSD
SDSal[0,4]=Box5SalSD
SDSal[0,5]=Box6SalSD
SDSal[0,6]=Box7SalSD
SDSal=np.round(SDSal,2)

SDPress=np.zeros([1,7])
SDPress[0,0]=Box1PressSD
SDPress[0,1]=Box2PressSD
SDPress[0,2]=Box3PressSD
SDPress[0,3]=Box4PressSD
SDPress[0,4]=Box5PressSD
SDPress[0,5]=Box6PressSD
SDPress[0,6]=Box7PressSD
SDPress=np.round(SDPress,2)

SDOx=np.zeros([1,7])
SDOx[0,0]=Box1OxSD
SDOx[0,1]=Box2OxSD
SDOx[0,2]=Box3OxSD
SDOx[0,3]=Box4OxSD
SDOx[0,4]=Box5OxSD
SDOx[0,5]=Box6OxSD
SDOx[0,6]=Box7OxSD
SDOx=np.round(SDOx,1)

SDSi=np.zeros([1,7])
SDSi[0,0]=Box1SiSD
SDSi[0,1]=Box2SiSD
SDSi[0,2]=Box3SiSD
SDSi[0,3]=Box4SiSD
SDSi[0,4]=Box5SiSD
SDSi[0,5]=Box6SiSD
SDSi[0,6]=Box7SiSD
SDSi=np.round(SDSi,1)

SDpCFC12=np.zeros([1,7])
SDpCFC12[0,0]=Box1pCFC12SD
SDpCFC12[0,1]=Box2pCFC12SD
SDpCFC12[0,2]=Box3pCFC12SD
SDpCFC12[0,3]=Box4pCFC12SD
SDpCFC12[0,4]=Box5pCFC12SD
SDpCFC12[0,5]=Box6pCFC12SD
SDpCFC12[0,6]=Box7pCFC12SD
SDpCFC12=np.round(SDpCFC12,6)

# Calculate carbonate ion concentration from the GLODAP data-------------------

# Setup
BT=1.179e-5*Sal
IonS = 19.924*Sal/(1000-1.005*Sal)

#Calculate carbonate
K0 = np.exp((A1_C + A2_C * (100.0/((Temp+Kelv))) + A3_C*np.log((Temp+Kelv)/100.0) + (Sal) * (B1_C + B2_C *((Temp+Kelv)/100.0) + B3_C * (((Temp+Kelv)/100.0)**2)))); # Weiss 1974
K1 = 10**-((3633.86/(Temp+Kelv)) - 61.2172 + (9.67770 * np.log(Temp+Kelv)) - (0.011555 * Sal) + (0.0001152 * (Sal**2))); #Lueker et al 2000
K2 = 10**-((471.78/((Temp+Kelv))) + 25.9290 - (3.16967 * np.log((Temp+Kelv))) - (0.01781 * Sal) + (0.0001122 * (Sal**2))); # Lueker et al 2000
KB = np.exp((-8966.9-2890.53*(Sal**0.5)-77.942*Sal+1.728*(Sal**1.5)-0.0996*(Sal**2))/(Temp+Kelv)+(148.0248+137.1942*(Sal**0.5)+(1.62142*Sal))+(-24.4344-25.085*(Sal**0.5)-0.2474*Sal)*np.log((Temp+Kelv))+0.053105*(Sal**0.5)*(Temp+Kelv)) #Dickson 1990
KW = np.exp(148.96502+(-13847.26/(Temp+Kelv))-(23.6521*np.log((Temp+Kelv)))+((Sal**0.5)*(-5.977+118.67/(Temp+Kelv)+1.0495*np.log((Temp+Kelv))))-0.01615*Sal) #Millero 1995
KS = np.exp(-8904.2/(Temp+Kelv)+117.4-19.334*np.log((Temp+Kelv))+(-458.79/(Temp+Kelv)+3.5913)*IonS**0.5+(188.74/(Temp+Kelv)-1.5998)*IonS+(-12.1652/(Temp+Kelv)+0.07871)*IonS*IonS)*(1-0.001005*Sal)
KP1 = np.exp(-4576.752/(Temp+Kelv)+115.54-18.453*np.log((Temp+Kelv))+(-106.736/(Temp+Kelv)+0.69171)*Sal**0.5+(-0.65643/(Temp+Kelv)-0.01844)*Sal) #Yao and Millero, 1995
KP2 = np.exp(-8814.715/(Temp+Kelv)+172.1033-27.927*np.log((Temp+Kelv))+(-160.34/(Temp+Kelv)+1.3566)*Sal**0.5+(0.37335/(Temp+Kelv)-0.05778)*Sal) #Yao and Millero, 1995
KP3 = np.exp(-3070.75/(Temp+Kelv)-18.126+(17.27039/(Temp+Kelv)+2.81197)*Sal**0.5+(-44.99486/(Temp+Kelv)-0.09984)*Sal) #Yao and Millero, 1995
KF = np.exp(1590.2/(Temp+Kelv)-12.641+1.525*IonS**0.5)*(1-0.001005*Sal)

# Create pH components

# calculations (Follows et al, 2006)

BAlk=BT*KB/(H+KB)
denom=H*H*H+(KP1*H*H)+(KP1*KP2*H)+(KP1*KP2*KP3)
SiAlk=(Si*1e-6)*KS/(KS+H)
H3PO4g=((Phos*1e-6)*H*H*H)/denom
H2PO4g=((Phos*1e-6)*KP1*H*H)/denom
HPO4g=((Phos*1e-6)*KP1*KP2*H)/denom
PO4g=((Phos*1e-6)*KP1*KP2*KP3)/denom
fg = -BAlk-(KW/H)+H-HPO4g-2*PO4g+H3PO4g-SiAlk

CAlk=Alk*1e-6+fg # convert umol/kg into mol/kg by multiplying by 1e-6
gamma=TCO2*1e-6/CAlk # convert umol/kg into mol/kg by multiplying by 1e-6
DM=(1-gamma)*(1-gamma)*K1*K1-4*K1*K2*(1-2*gamma)
H=0.5*((gamma-1)*K1+(DM)**0.5)
pCO2 = (TCO2*1e-6/K0)*((H**2)/((H**2)+K1*H+K1*K2))*1e6 # convert umol/kg into mol/kg by multiplying by 1e-6
pH = np.log(H)/np.log(0.1)
H2CO3 = pCO2*K0
HCO3 = (H2CO3*K1)/H
CO23 = (HCO3*K2)/H
DIC = H2CO3+HCO3+CO23
CO23=np.round(CO23,0)

## Print output - comment out if not using

#print(Phos)
#print(TCO2)
#print(Alk)
#print(d14C)
#print(d13C)
#print(Temp)
#print(Sal)
#print(Press)
#print(Ox)
#print(Si)
#print(pCFC12)
#print(CO23)


##Save text files (optional) - comment out if not using

#np.savetxt(DataDir+'GLODAP_Phos.txt',Phos)
#np.savetxt(DataDir+'GLODAP_DIC.txt',TCO2)
#np.savetxt(DataDir+'GLODAP_Alk.txt',Alk)
#np.savetxt(DataDir+'GLODAP_d14C.txt',d14C)
#np.savetxt(DataDir+'GLODAP_d13C.txt',d13C)
#np.savetxt(DataDir+'GLODAP_Temp.txt',Temp)
#np.savetxt(DataDir+'GLODAP_Sal.txt',Sal)
#np.savetxt(DataDir+'GLODAP_Press.txt',Press)
#np.savetxt(DataDir+'GLODAP_Ox.txt',Ox)
#np.savetxt(DataDir+'GLODAP_Si.txt',Si)
#np.savetxt(DataDir+'GLODAP_pCFC12.txt',pCFC12)
#np.savetxt(DataDir+'GLODAP_CO23.txt',CO23)
#
#np.savetxt(DataDir+'GLODAP_PhosSD.txt',SDPhos)
#np.savetxt(DataDir+'GLODAP_DICSD.txt',SDTCO2)
#np.savetxt(DataDir+'GLODAP_AlkSD.txt',SDAlk)
#np.savetxt(DataDir+'GLODAP_d14CSD.txt',SDd14C)
#np.savetxt(DataDir+'GLODAP_d13CSD.txt',SDd13C)
#np.savetxt(DataDir+'GLODAP_TempSD.txt',SDTemp)
#np.savetxt(DataDir+'GLODAP_SalSD.txt',SDSal)
#np.savetxt(DataDir+'GLODAP_PressSD.txt',SDPress)
#np.savetxt(DataDir+'GLODAP_OxSD.txt',SDOx)
#np.savetxt(DataDir+'GLODAP_SiSD.txt',SDSi)
#np.savetxt(DataDir+'GLODAP_pCFC12SD.txt',SDpCFC12)




