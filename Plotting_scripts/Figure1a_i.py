#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure 1 a-i

Necessary datasets (ESGF):
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.epc100.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.fgco2.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.pp.gn

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import cmocean as cm
import xesmf as xe

#%% Gas exchange plots (Figure 1a, b and c)
data1='/Users/daan/CESM2_data'   # Location of dataset(s) 
var1='fgco2'                     # Variable name

#%% Call on datasets (might need to change name of dataset)
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Load data
l1=load_var1[var1][:,:,:].compute().squeeze()
l2=load_var2[var1][:,:,:].compute().squeeze()

#%% Concatenate datasets in time
VAR=xr.concat([l1,l2],dim='time')

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
VAR_gr=np.roll(VAR_gr,-180)

#%% Load longitude and latitude data for 'gr' grid
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Plotting variables
FS=20

vmn1=-1e-9
vmx1=1e-9

vmn2=-1e-9
vmx2=1e-9

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='BrBG_r')
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Gas exchange [kg C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = False
gl.left_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='BrBG_r')
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Gas exchange [kg C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0)-np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$Gas exchange [kg C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

#plt.savefig('figure1_1.png', format='png', dpi=300)

#%% Net primary production plots (Figure 1 d, e and f)
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='pp'                        # Variable name

#%% Call on datasets
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Load data and integrate over top 100m
l1=load_var1[var1][:,:10,:,:].sum('lev').compute().squeeze()*10
l2=load_var2[var1][:,:10,:,:].sum('lev').compute().squeeze()*10

#%% Concatenate datasets in time
VAR=xr.concat([l1,l2],dim='time')

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
VAR_gr=np.roll(VAR_gr,-180)

#%% Load longitude and latitude data for 'gr' grid
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Plotting variables
FS=20

vmn1=0
vmx1=1e-6

vmn2=-1e-7
vmx2=1e-7

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0)-np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$NPP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Save figure as...
#plt.savefig('figure1_2.png', format='png', dpi=300)

#%% Export production plots (Figure 1 g, h and i)
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='epc100'                    # Variable name

#%% Call on dataset
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_2.nc')

#%% Load data
l1=load_var1[var1][:,:,:].compute().squeeze()
l2=load_var2[var1][:,:,:].compute().squeeze()

#%% Concatenate datasets in time
VAR=xr.concat([l1,l2],dim='time')

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
VAR_gr=np.roll(VAR_gr,-180)

#%% Load longitude and latitude data for 'gr' grid
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Plotting variables
FS=20               # Fontsize 

vmn1=0              # Bounds for colormap
vmx1=1e-7

vmn2=-5e-8          # Bounds for colormap
vmx2=5e-8

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.matter)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('EP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.matter)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('EP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0)-np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$EP [mol C m$^{-2}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
gl.top_labels = False
gl.right_labels = False
gl.rotate_labels = False

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

# Save figure as...
#plt.savefig('figure1_3.png', format='png', dpi=300)
