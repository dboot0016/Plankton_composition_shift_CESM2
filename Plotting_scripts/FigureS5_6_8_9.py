"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Supplementary plots
Script for plotting Figure S5, S6, S8, S9

Necessary datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.thetao.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.so.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.mlotstmax.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.pp.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.epc100.gn

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
import gsw
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#%% Sea surface temperature
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='thetao'    # Variable name

#%% Call on data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Load data and select surface layer
l1=load_var1[var1][:,0,:,:].compute().squeeze()
l2=load_var2[var1][:,0,:,:].compute().squeeze()

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
vmx1=35

vmn2=-5
vmx2=5

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.thermal)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('SST [$^{\circ}$C]', fontsize=16)
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
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.thermal)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('SST [$^{\circ}$C]', fontsize=16)
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
cbar.ax.set_xlabel('$\Delta$SST [$^{\circ}$C]', fontsize=16)
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
#plt.savefig('figureS5.png', format='png', dpi=300)

#%% Maximum mixed layer depth
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='mlotstmax'                 # Variable name

#%% Call on data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_2.nc')

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

vmn1=0
vmx1=300

vmn2=-200
vmx2=200

#%% Figure 
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.deep_r)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('MLD (max) [m]', fontsize=16)
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
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.deep_r)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('MLD (max) [m]', fontsize=16)
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
cbar.ax.set_xlabel('$\Delta$MLD (max) [m]', fontsize=16)
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
#plt.savefig('figureS6.png', format='png', dpi=300)

#%% Stratification plots (Figure S6a, b and c)
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='thetao'                    # Variable name
var2='so'

#%% Call on datasets
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

load_var3 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_1.nc')
load_var4 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_2.nc')

#%% Temperature at 0m and 200m
l1=load_var1[var1][:,0,:,:].compute().squeeze()
l2=load_var2[var1][:,0,:,:].compute().squeeze()

l3=load_var1[var1][:,20,:,:].compute().squeeze()
l4=load_var2[var1][:,20,:,:].compute().squeeze()

#%% Salinity at 0m and 200m
l5=load_var3[var2][:,0,:,:].compute().squeeze()
l6=load_var4[var2][:,0,:,:].compute().squeeze()

l7=load_var3[var2][:,20,:,:].compute().squeeze()
l8=load_var4[var2][:,20,:,:].compute().squeeze()

#%% Concatenate datasets in time
VAR1=xr.concat([l1,l2],dim='time')
VAR2=xr.concat([l3,l4],dim='time')
VAR3=xr.concat([l5,l6],dim='time')
VAR4=xr.concat([l7,l8],dim='time')

#%% Determine densities at 0m and 200m and difference between the two
CT1=gsw.CT_from_pt(VAR3,VAR1)
p1=gsw.p_from_z(0,VAR3['lat'])

rho1=gsw.density.rho(VAR3,CT1,p1)

CT2=gsw.CT_from_pt(VAR4,VAR2)
p2=gsw.p_from_z(-VAR4['lev']/100,VAR4['lat'])

rho2=gsw.density.rho(VAR4,CT2,p2)

VAR=rho2-rho1

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
vmx1=10

vmn2=-5
vmx2=5

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15*12,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.dense_r)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Stratification [kg m$^{-3}$]', fontsize=16)
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
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15*12:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.dense_r)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Stratification [kg m$^{-3}$]', fontsize=16)
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
cbar.ax.set_xlabel('$\Delta$Stratification [kg m$^{-3}$]', fontsize=16)
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
#plt.savefig('figureS8.png', format='png', dpi=300)

#%% E-ratio: call on NPP
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='pp'    # Variable name

#%% Call on NPP data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Load in data and integrate over top 100m
l1=load_var1[var1][:,:10,:].sum('lev').compute().squeeze()*10
l2=load_var2[var1][:,:10,:].sum('lev').compute().squeeze()*10

#%% Concatenate datasets in time
VAR=xr.concat([l1,l2],dim='time')

#%% E-ratio: call on EP data
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='epc100'    # Variable name

#%% Call on EP data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_2.nc')

#%% Load in datasets
l1=load_var1[var1][:,:,:].compute().squeeze()
l2=load_var2[var1][:,:,:].compute().squeeze()

#%% Concatenate datasets in time
VAR2=xr.concat([l1,l2],dim='time')

#%% Take yearly averages of EP data
varr2=np.zeros((86,384,320))
for i in range(86):
    varr2[i,:,:]=VAR2[i*12:(i+1)*12,:,:].mean('time')

#%% Determine e-ratio
VAR2=varr2*xr.ones_like(VAR)
VARe=VAR2/VAR

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VARe, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VARe)
VAR_gr=np.roll(VAR_gr,-180)

#%% Load longitude and latitude data for 'gr' grid
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%% Plotting variables
FS=20

vmn1=0.075
vmx1=0.275

vmn2=-0.075
vmx2=0.075

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
ax.set_extent([-100,20,10,80])
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdYlBu_r')
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('e-ratio [-]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
ax.set_extent([-100,20,10,80])
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap='RdYlBu_r')
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('e-ratio [-]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
ax.set_extent([-100,20,10,80])
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0)-np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

# Colorbar specifics
cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$e-ratio [-]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Grid lines and longitude and latitude notations
gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

# Save figure as...
#plt.savefig('figureS9.png', format='png', dpi=300)