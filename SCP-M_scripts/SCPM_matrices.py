#!/usr/bin/python3 # 

# [simple carbon project] model SCP-M
# 7 ocean box model calculations (Box 0 - Box 6) plus atmosphere, zonally averaged
# The Australian National University

# O'Neill, C.M. and A. McC. Hogg, M. Ellwood, B.N. Opydke and S.M. Eggins
# The [simple carbon project] model, in process


## SCP-M matrices file
# This module controls the matrices that govern the operation of key flux parameters 
# in the SCP-M model: ocean physical, biological, carbonate rain and piston velocity
# Disclaimer: SCP-M comes with no warranty or guarantee, use at your own risk etc (have to say this apparently)

#Import modules
import math
import numpy as np
import matplotlib.pyplot as plt

from SCPM_parameters import (Varr,Fract,Zed,bScale,PV0,PV1,PV4,PV6,gamma2,secsday, 
F0,F1,F4,F6,Psi1,Psi2,gamma1,rcp,rcpFe,rcpSi,rcpO,FCA)

# Ocean physical---------------------------------------------------------------

def makeT1mat(Psi1, Varr, Fract):
    makeT1mat = np.zeros([7,7])
    makeT1mat[3,3]=-1
    makeT1mat[3,5]=1
    makeT1mat[4,4]=-1
    makeT1mat[4,6]=Fract
    makeT1mat[4,3]=1-Fract
    makeT1mat[5,4]=1
    makeT1mat[5,5]=-1
    makeT1mat[6,3]=Fract
    makeT1mat[6,6]=-Fract
    return makeT1mat
    
def makeT2mat(Psi2, Varr):
    makeT2mat = np.zeros([7,7])      
    makeT2mat[1,1]=-1
    makeT2mat[1,2]=1
    makeT2mat[2,2]=-1
    makeT2mat[2,6]=1
    makeT2mat[3,1]=1
    makeT2mat[3,3]=-1
    makeT2mat[6,3]=1  
    makeT2mat[6,6]=-1
    return makeT2mat
    
def makeE1mat(gamma1,Varr):
    makeE1mat = np.zeros([7,7])
    makeE1mat[3,3]=-1
    makeE1mat[3,5]=1    
    makeE1mat[5,3]=1
    makeE1mat[5,5]=-1       
    return makeE1mat

def makeE2mat(gamma2,Varr):
    makeE2mat = np.zeros([7,7])
    makeE2mat[0,0]=-1
    makeE2mat[0,2]=1
    makeE2mat[2,0]=1
    makeE2mat[2,2]=-1
    return makeE2mat

# Marine biological flux-------------------------------------------------------

def makeBiomat_out(Zed, Varr):
    makeBiomat_out=np.zeros([7,7])
    makeBiomat_out[0,0]=-1
    makeBiomat_out[1,1]=-1
    makeBiomat_out[2,0]=-1
    makeBiomat_out[3,0]=-1
    makeBiomat_out[3,1]=-1
    makeBiomat_out[3,6]=-1
    makeBiomat_out[4,4]=-1
    makeBiomat_out[5,0]=-0 # assuming zero soft biological flux from the abyssal
    makeBiomat_out[5,1]=-0 # ocean
    makeBiomat_out[5,4]=-0
    makeBiomat_out[6,6]=-1
    return makeBiomat_out

def makeBiomat_in(Zed, Varr):
    makeBiomat_in=np.zeros([7,7])
    makeBiomat_in[2,0]=1
    makeBiomat_in[3,0]=1
    makeBiomat_in[3,1]=1
    makeBiomat_in[3,6]=1
    makeBiomat_in[5,0]=1
    makeBiomat_in[5,1]=1
    makeBiomat_in[5,4]=1
    makeBiomat_in[5,6]=1
    return makeBiomat_in

#  Martin surface box biological carbon flux values at 100m depth in mols C m-2 yr-1

def makeZarr(Zed):
    makeZarr=np.zeros([7,1]) 
    makeZarr[0,0]=Zed
    makeZarr[1,0]=Zed
    makeZarr[4,0]=Zed
    makeZarr[6,0]=Zed
    return makeZarr


# Air to sea surface gas transfer----------------------------------------------
def makeFXarr(PV0, PV1,PV4,PV6, secsday):
    makeFXarr = np.zeros([7,1])
    makeFXarr[0,0] = (PV0/secsday)
    makeFXarr[1,0] = (PV1/secsday)
    makeFXarr[4,0] = (PV4/secsday)
    makeFXarr[6,0] = (PV6/secsday)
    return makeFXarr










