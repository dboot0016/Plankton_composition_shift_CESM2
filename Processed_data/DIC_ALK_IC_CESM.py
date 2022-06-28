#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 09:08:56 2022

@author: daan
"""
import xarray as xr 
import numpy as np
import xesmf as xe
import matplotlib.pyplot as plt

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='dissic'    # Variable name
var2='talk'

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1y_1.nc')

#%%
l1=load_var1[var1][0,:,:,:].squeeze()
l2=load_var2[var2][0,:,:,:].squeeze()

#%%
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(l1, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(l1)
L1=np.roll(VAR_gr,0)

regridder = xe.Regridder(l2, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(l2)
L2=np.roll(VAR_gr,0)

#%%
D=l1['lev']/100

dz = (l1.lev.shift(lev=-1)-l1.lev.shift(lev=1))/2
dz[0]  = l1.lev[0] + (l1.lev[1]-l1.lev[0])/2
dz[-1] = 1.5*(l1.lev[-1]-l1.lev[-2])

dz=dz/100

#%% Load in area
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%%
y=(np.isnan(L1)==0)

#%%
Vol=np.zeros((60,180,360))
for i in range(180):
    for j in range(360):
        Vol[:,i,j]=y[:,i,j]*np.array(area)[i,j]*dz
        
#%%
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#%%
Box1_d=100
Box1_1=(find_nearest(D, Box1_d))
Box1_2=-40+90
Box1_3=40+90

Box2_d=250
Box2_1=(find_nearest(D, Box2_d))
Box2_2=40+90
Box2_3=60+90

Box3_d=1000
Box3_1=(find_nearest(D, Box3_d))
Box3_2=-40+90
Box3_3=40+90

Box5_d=2500
Box5_1=(find_nearest(D, Box5_d))
Box5_2=-75+90
Box5_3=-40+90

Box6_d=4000
Box6_1=(find_nearest(D, Box6_d))
Box6_2=-75+90
Box6_3=60+90

Box7_d=250
Box7_1=(find_nearest(D, Box7_d))
Box7_2=-60+90
Box7_3=-40+90

#%%
DIC1=np.nansum(L1[:Box1_1,Box1_2:Box1_3,:]*Vol[:Box1_1,Box1_2:Box1_3,:])/np.nansum(Vol[:Box1_1,Box1_2:Box1_3,:])
DIC2=np.nansum(L1[:Box2_1,Box2_2:Box2_3,:]*Vol[:Box2_1,Box2_2:Box2_3,:])/np.nansum(Vol[:Box2_1,Box2_2:Box2_3,:])
DIC3=np.nansum(L1[Box1_1:Box3_1,Box3_2:Box3_3,:]*Vol[Box1_1:Box3_1,Box3_2:Box3_3,:])/np.nansum(Vol[Box1_1:Box3_1,Box3_2:Box3_3,:])
DIC5=np.nansum(L1[:Box5_1,Box5_2:Box5_3,:]*Vol[:Box5_1,Box5_2:Box5_3,:])/np.nansum(Vol[:Box5_1,Box5_2:Box5_3,:])
DIC6=np.nansum(L1[Box5_1:Box6_1,Box6_2:Box6_3,:]*Vol[Box5_1:Box6_1,Box6_2:Box6_3,:])/np.nansum(Vol[Box5_1:Box6_1,Box6_2:Box6_3,:])
DIC7=np.nansum(L1[:Box7_1,Box7_2:Box7_3,:]*Vol[:Box7_1,Box7_2:Box7_3,:])/np.nansum(Vol[:Box7_1,Box7_2:Box7_3,:])

ALK1=np.nansum(L2[:Box1_1,Box1_2:Box1_3,:]*Vol[:Box1_1,Box1_2:Box1_3,:])/np.nansum(Vol[:Box1_1,Box1_2:Box1_3,:])
ALK2=np.nansum(L2[:Box2_1,Box2_2:Box2_3,:]*Vol[:Box2_1,Box2_2:Box2_3,:])/np.nansum(Vol[:Box2_1,Box2_2:Box2_3,:])
ALK3=np.nansum(L2[Box1_1:Box3_1,Box3_2:Box3_3,:]*Vol[Box1_1:Box3_1,Box3_2:Box3_3,:])/np.nansum(Vol[Box1_1:Box3_1,Box3_2:Box3_3,:])
ALK5=np.nansum(L2[:Box5_1,Box5_2:Box5_3,:]*Vol[:Box5_1,Box5_2:Box5_3,:])/np.nansum(Vol[:Box5_1,Box5_2:Box5_3,:])
ALK6=np.nansum(L2[Box5_1:Box6_1,Box6_2:Box6_3,:]*Vol[Box5_1:Box6_1,Box6_2:Box6_3,:])/np.nansum(Vol[Box5_1:Box6_1,Box6_2:Box6_3,:])
ALK7=np.nansum(L2[:Box7_1,Box7_2:Box7_3,:]*Vol[:Box7_1,Box7_2:Box7_3,:])/np.nansum(Vol[:Box7_1,Box7_2:Box7_3,:])

#%%
DIC4_1=np.nansum(L1[Box3_1:Box5_1,Box7_2:Box2_3,:]*Vol[Box3_1:Box5_1,Box7_2:Box2_3,:])
DIC4_2=np.nansum(L1[Box2_1:Box3_1,Box7_2:Box7_3,:]*Vol[Box2_1:Box3_1,Box7_2:Box7_3,:])
DIC4_3=np.nansum(L1[Box2_1:Box3_1,Box2_2:Box2_3,:]*Vol[Box2_1:Box3_1,Box2_2:Box2_3,:])
V4_3=np.nansum(Vol[Box3_1:Box5_1,Box7_2:Box2_3,:])+np.nansum(Vol[Box2_1:Box3_1,Box7_2:Box7_3,:])+np.nansum(Vol[Box2_1:Box3_1,Box2_2:Box2_3,:])

DIC4=(DIC4_1+DIC4_2+DIC4_3)/V4_3

ALK4_1=np.nansum(L2[Box3_1:Box5_1,Box7_2:Box2_3,:]*Vol[Box3_1:Box5_1,Box7_2:Box2_3,:])
ALK4_2=np.nansum(L2[Box2_1:Box3_1,Box7_2:Box7_3,:]*Vol[Box2_1:Box3_1,Box7_2:Box7_3,:])
ALK4_3=np.nansum(L2[Box2_1:Box3_1,Box2_2:Box2_3,:]*Vol[Box2_1:Box3_1,Box2_2:Box2_3,:])

ALK4=(ALK4_1+ALK4_2+ALK4_3)/V4_3

#%%
DIC=np.array([DIC1,DIC2,DIC3,DIC4,DIC5,DIC6,DIC7])
ALK=np.array([ALK1,ALK2,ALK3,ALK4,ALK5,ALK6,ALK7])