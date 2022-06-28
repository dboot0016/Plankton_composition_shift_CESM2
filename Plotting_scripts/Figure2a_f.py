#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, D., Von Der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure 2 a-f

Necessary datasets (ESGF):
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limirrdiat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limirrpico.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limndiat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limnpico.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limfediat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limfediat.gn
        
Necessary datasets (processed):
- limsidiat_esm585_1.nc
- limpdiat_esm585_1.nc
- limppico_esm585_1.nc
- totD_esm585_100_1.nc
- totP_esm585_100_1.nc
- thetao_esm585_zint_1.nc

@author: Daan Boot (d.boot@uu.nl)
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

#%% NPP plots per PFT (Figure 2a-f)
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='limnpico'                  # Nitrogen limitation small phytoplankton
var2='limndiat'                  # Nitrogen limitation diatoms
var3='limfepico'                 # Iron limitation small phytoplankton
var4='limfediat'                 # Iron limitation diatoms
var5='limirrpico'                # Light limitation small phytoplankton
var6='limirrdiat'                # Light limitation diatoms
var7='thetao'                    # Potential temperature
var8='limsidiat'                 # Silicate limitation diatoms
var9='limpdiat'                  # Phosphorus limitation diatoms
var10='limppico'                 # Phosphorus limitation small phytoplankton
var11='totD'                     # Integrated diatom biomass (top 100m)
var12='totP'                     # Integrated small phytoplankton biomass (top 100m)
            
#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')
load_var3 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_1.nc')
load_var4 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_2.nc')
load_var5 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_1.nc')
load_var6 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_2.nc')

load_var7 = xr.open_dataset(f'{data1}/'+var4+'_esm585_1_1.nc')
load_var8 = xr.open_dataset(f'{data1}/'+var4+'_esm585_1_2.nc')
load_var9 = xr.open_dataset(f'{data1}/'+var5+'_esm585_1_1.nc')
load_var10 = xr.open_dataset(f'{data1}/'+var5+'_esm585_1_2.nc')
load_var11 = xr.open_dataset(f'{data1}/'+var6+'_esm585_1_1.nc')
load_var12 = xr.open_dataset(f'{data1}/'+var6+'_esm585_1_2.nc')

#%%
load_var13 = xr.open_dataset(f'{data1}/'+var7+'_esm585_zint_1.nc')
load_var14 = xr.open_dataset(f'{data1}/'+var8+'_esm585_1.nc')
load_var15 = xr.open_dataset(f'{data1}/'+var9+'_esm585_1.nc')
load_var16 = xr.open_dataset(f'{data1}/'+var10+'_esm585_1.nc')
load_var17 = xr.open_dataset(f'{data1}/'+var11+'_esm585_100_1.nc')
load_var18 = xr.open_dataset(f'{data1}/'+var12+'_esm585_100_1.nc')

#%%
l1=load_var1[var1].compute().squeeze()
l2=load_var2[var1].compute().squeeze()
l3=load_var3[var2].compute().squeeze()
l4=load_var4[var2].compute().squeeze()
l5=load_var5[var3].compute().squeeze()
l6=load_var6[var3].compute().squeeze()
l7=load_var7[var4].compute().squeeze()
l8=load_var8[var4].compute().squeeze()
l9=load_var9[var5].compute().squeeze()
l10=load_var10[var5].compute().squeeze()
l11=load_var11[var6].compute().squeeze()
l12=load_var12[var6].compute().squeeze()

#%%
varx='__xarray_dataarray_variable__'
l13=load_var13[var7].compute().squeeze()
l14=load_var14[varx].compute().squeeze()
l15=load_var15[varx].compute().squeeze()
l16=load_var16[varx].compute().squeeze()
l17=load_var17['po4'].compute().squeeze()
l18=load_var18['po4'].compute().squeeze()

#%%
limN_pico=xr.concat([l1,l2],dim='time')
limN_diat=xr.concat([l3,l4],dim='time')
limF_pico=xr.concat([l5,l6],dim='time')
limF_diat=xr.concat([l7,l8],dim='time')
limI_pico=xr.concat([l9,l10],dim='time')
limI_diat=xr.concat([l11,l12],dim='time')
limP_pico=l16
limP_diat=l15
limS_diat=l14

theta=l13
si=l14/150
p_diat=l17/100
p_pico=l18/100

#%%
Tf=1.7**((theta-30)/10)
uref_p=5/(86400)
uref_di=5/(86400)

#%%
limf_P=np.array(limF_pico)
limn_P=np.array(limN_pico)
limp_p=np.array(limP_pico)

limf_Di=np.array(limF_diat)
limn_Di=np.array(limN_diat)
lims_di=np.array(limS_diat)
limp_di=np.array(limP_diat)

limI_Pico=np.array(limI_pico)
limI_Diat=np.array(limI_diat)

tf=np.array(Tf)

#%%
limf_p=np.zeros((86,384,320))
limn_p=np.zeros((86,384,320))

limf_di=np.zeros((86,384,320))
limn_di=np.zeros((86,384,320))

limf_dz=np.zeros((86,384,320))
limt_dz=np.zeros((86,384,320))

limi_pico=np.zeros((86,384,320))
limi_diat=np.zeros((86,384,320))
limi_diaz=np.zeros((86,384,320))

TF=np.zeros((86,384,320))

for i in range(86):
    limf_p[i,:,:]=np.mean(limf_P[i*12:(i+1)*12,:,:],axis=0)
    limn_p[i,:,:]=np.mean(limf_P[i*12:(i+1)*12,:,:],axis=0)
    
    limf_di[i,:,:]=np.mean(limf_Di[i*12:(i+1)*12,:,:],axis=0)
    limn_di[i,:,:]=np.mean(limn_Di[i*12:(i+1)*12,:,:],axis=0)
    
    limi_pico[i,:,:]=np.mean(limI_Pico[i*12:(i+1)*12,:,:],axis=0)
    limi_diat[i,:,:]=np.mean(limI_Diat[i*12:(i+1)*12,:,:],axis=0)
    
    TF[i,:,:]=np.mean(tf[i*12:(i+1)*12,:,:],axis=0)
    
#%%
limp=np.zeros((86,384,320))
limdi=np.zeros((86,384,320))
limdz=np.zeros((86,384,320))

for i in range(384):
    for j in range(320):
        for k in range(86):
            limp[k,i,j]=min([limn_p[k,i,j],limp_p[k,i,j],limf_p[k,i,j]])
            limdi[k,i,j]=min([limn_di[k,i,j],limp_di[k,i,j],limf_di[k,i,j],lims_di[k,i,j]])
               
#%%
for i in range (384):
    for j in range(320):
        if np.isnan(limn_p[0,i,j]):
            limp[:,i,j]='nan'
            limdi[:,i,j]='nan'

#%%
npp_p=uref_p*TF*limi_pico*limp*p_pico
npp_di=uref_di*TF*limi_diat*limdi*p_diat

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var9+'_esm585_1.nc')

#%%
loaded_var=load_var1[varx][:,:,:].compute().squeeze()

#%%
npp_p=npp_p*xr.ones_like(loaded_var)
npp_di=npp_di*xr.ones_like(loaded_var)

#%%
npp_p1=npp_p[:15,:,:].mean('time')
npp_p2=npp_p[-15:,:,:].mean('time')
dnpp_p=npp_p2-npp_p1

npp_di1=npp_di[:15,:,:].mean('time')
npp_di2=npp_di[-15:,:,:].mean('time')
dnpp_di=npp_di2-npp_di1

#%%
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(npp_p1, ds_out, 'bilinear',periodic=True)
NPP_p1 = regridder(npp_p1)
NPP_p1=np.roll(NPP_p1,-180)

regridder = xe.Regridder(npp_p2, ds_out, 'bilinear',periodic=True)
NPP_p2 = regridder(npp_p2)
NPP_p2=np.roll(NPP_p2,-180)

regridder = xe.Regridder(dnpp_p, ds_out, 'bilinear',periodic=True)
dNPP_p = regridder(dnpp_p)
dNPP_p=np.roll(dNPP_p,-180)

regridder = xe.Regridder(npp_di1, ds_out, 'bilinear',periodic=True)
NPP_di1 = regridder(npp_di1)
NPP_di1=np.roll(NPP_di1,-180)

regridder = xe.Regridder(npp_di2, ds_out, 'bilinear',periodic=True)
NPP_di2 = regridder(npp_di2)
NPP_di2=np.roll(NPP_di2,-180)

regridder = xe.Regridder(dnpp_di, ds_out, 'bilinear',periodic=True)
dNPP_di = regridder(dnpp_di)
dNPP_di=np.roll(dNPP_di,-180)

#%%
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()

#%%
FS=20

vmn1=0
vmx1=7.5e-10

vmn2=-7.5e-10
vmx2=7.5e-10

vmn3=0
vmx3=1e-11

vmn4=-1e-11
vmx4=1e-11

level1=np.linspace(vmn1,vmx1,7)
level2=np.linspace(vmn2,vmx2,7)

fig = plt.figure(figsize=(20, 5))

ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,NPP_di1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (di) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

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

ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,NPP_di2,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (di) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

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

ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,dNPP_di,vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (di) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

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

#plt.savefig('figure2_1.png', format='png', dpi=300)

#%%
fig = plt.figure(figsize=(20, 5))

ax = fig.add_subplot(1, 3, 1, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,NPP_p1,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2015)+ ' - ' +str(2030),fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (sp) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

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

ax = fig.add_subplot(1, 3, 2, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,NPP_p2,vmin=vmn1,vmax=vmx1,transform=ccrs.PlateCarree(),cmap=cm.cm.algae)
ax.set_title('Average ' + str(2086)+ ' - ' +str(2101),fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (sp) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

ax.text(-16000000,-8500000,'80$^{\circ}$S',fontsize=14)
ax.text(-20000000,-4800000,'40$^{\circ}$S',fontsize=14)
ax.text(-19000000,-500000,'0$^{\circ}$',fontsize=14)
ax.text(-20000000,4000000,'40$^{\circ}$N',fontsize=14)
ax.text(-16000000,7700000,'80$^{\circ}$N',fontsize=14)

ax.text(-10000000,-10000000,'180$^{\circ}$E',fontsize=14)
ax.text(-4000000,-10000000,'270$^{\circ}$E',fontsize=14)
ax.text(2000000,-10000000,'0$^{\circ}$E',fontsize=14)
ax.text(7000000,-10000000,'90$^{\circ}$E',fontsize=14)

ax = fig.add_subplot(1, 3, 3, projection=ccrs.Robinson(-60))
ax.coastlines(resolution='50m',zorder=5)
ax.add_feature(cfeature.LAND,zorder=4)
im=plt.pcolormesh(lon,lat,dNPP_p,vmin=vmn2,vmax=vmx2,transform=ccrs.PlateCarree(),cmap='RdBu_r')
ax.set_title('Difference',fontsize=FS)

cbar=plt.colorbar(im,orientation='horizontal',) 
cbar.ax.set_xlabel('NPP (sp) [mol C m$^{-3}$ s$^{-1}$]', fontsize=16)
cbar.ax.set_yticks(fontsize=16)
cbar.ax.tick_params(labelsize=16)
cbar.ax.xaxis.offsetText.set_fontsize(16)

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

#plt.savefig('figure2_2.png', format='png', dpi=300)