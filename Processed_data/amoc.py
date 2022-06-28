"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S15a-f

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.thetao.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.so.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.msftmz.gn

@author: Amber Boot (d.boot@uu.nl)
"""

import intake
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import cmocean as cm

import regionmask
from cmip6_preprocessing.regionmask import merged_mask
from cftime import DatetimeNoLeap
import xesmf as xe
import gsw

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error, r2_score

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='thetao'    # Variable name
var2='so'    # Variable name
var3='msftmz'    # Variable na

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%%
T1=load_var1[var1][:,:,271,:25].compute().squeeze()
T2=load_var1[var1][:,:,271,-36:].compute().squeeze()

TA=xr.concat([T2,T1],dim='nlon')

#%%
T1=load_var2[var1][:,:,271,:25].compute().squeeze()
T2=load_var2[var1][:,:,271,-36:].compute().squeeze()

TB=xr.concat([T2,T1],dim='nlon')

#%%
T=xr.concat([TA.mean('nlon'),TB.mean('nlon')],dim='time')

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_2.nc')

#%%
S1=load_var1[var2][:,:,271,:25].compute().squeeze()
S2=load_var1[var2][:,:,271,-36:].compute().squeeze()

SA=xr.concat([S2,S1],dim='nlon')

#%%
S1=load_var2[var2][:,:,271,:25].compute().squeeze()
S2=load_var2[var2][:,:,271,-36:].compute().squeeze()

SB=xr.concat([S2,S1],dim='nlon')

#%%
S=xr.concat([SA.mean('nlon'),SB.mean('nlon')],dim='time')

#%%
CT = gsw.conversions.CT_from_pt(S, T)
p = gsw.p_from_z(-T.lev/100,0)
rho = gsw.density.rho(S, CT,p)

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_2.nc')

#%%
l1=load_var1[var3][:,0,:,274].squeeze()
l2=load_var2[var3][:,0,:,274].squeeze()

M=xr.concat([l1,l2],dim='time')

#%%
M2=np.zeros((1032,60))
for i in range(60):
    M2[:,i]=M[:,i:i+1].mean('lev')

#%%
r=np.array(rho)
amoc=M2/r
AMOC=amoc*xr.ones_like(rho)

#%%
#AMOC.to_netcdf('amoc_esm585_1.nc')
