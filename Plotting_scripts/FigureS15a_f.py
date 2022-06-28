"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S15a-f

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.phydiat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.phypico.gn

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import cmocean as cm
import xesmf as xe
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#%% Load in biomass diatom
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='phydiat'                    # Variable name

#%% Call on dataset
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_2.nc')

#%% Load data + sum up over vertical
l1=load_var1[var1].sum(['lev']).compute().squeeze()*10
l2=load_var2[var1].sum(['lev']).compute().squeeze()*10

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
vmx1=1e-1

vmn2=-1e-1
vmx2=1e-1

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.haline)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Biomass (di) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.haline)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Biomass (di) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0)-np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$Biomass (di) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

#plt.savefig('figureS15d_f.png', format='png', dpi=300)

#%% Small phytoplankton biomass
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='phypico'                    # Variable name

#%% Call on dataset
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Load data + sum up over vertical
l1=load_var1[var1].sum(['lev']).compute().squeeze()*10
l2=load_var2[var1].sum(['lev']).compute().squeeze()*10

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
vmx1=0.5e-1

vmn2=-5e-2
vmx2=5e-2

#%% Figure
fig = plt.figure(figsize=(20, 5))

# Subplot 1
ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.haline)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Biomass (sp) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Subplot 2
ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0),vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.haline)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('Biomass (sp) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

# Subplot 3
ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,np.mean(VAR_gr[-15:,:,:],axis=0)-np.mean(VAR_gr[:15,:,:],axis=0),vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)
ax.set_extent([-100,20,10,80])

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False

gl.xlocator = mticker.FixedLocator([-100,-75,-50,-25, 0])
gl.ylocator = mticker.FixedLocator([0, 20, 40, 60,80])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('$\Delta$Biomass (sp) [mol C m$^{-2}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

#plt.savefig('figureS15a_c.png', format='png', dpi=300)

