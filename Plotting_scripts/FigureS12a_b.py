"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S14a-f

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
var9='limpdiat'                  # Nitrogen limitation small phytoplankton
var10='limppico'                 # Nitrogen limitation small phytoplankton
var11='totD'                     # Nitrogen limitation small phytoplankton
var12='totP'                     # Nitrogen limitation small phytoplankton
            
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

#%% Determine growth rate
mu_p=uref_p*TF*limi_pico*limp
mu_di=uref_di*TF*limi_diat*limdi

#%% Get a dataset with annual data (doesn't matter which)
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='phydiat'    # Variable name

load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1y_2.nc')

l1=load_var1[var1][:,0,:,:].compute().squeeze()
l2=load_var2[var1][:,0,:,:].compute().squeeze()

loaded_var=(xr.concat([l1,l2],dim='time'))

#%% Numpy array back to xarray
mu_p=mu_p*xr.ones_like(loaded_var)
mu_di=mu_di*xr.ones_like(loaded_var)

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(mu_p, ds_out, 'bilinear',periodic=True)
MU_P = regridder(mu_p)
MU_P=np.roll(MU_P,-180)

ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(mu_di, ds_out, 'bilinear',periodic=True)
MU_D = regridder(mu_di)
MU_D=np.roll(MU_D,-180)

#%% Load in area
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%% Multiply with area
mu_p=np.zeros((np.size(MU_P[:,0,0]),180,360))
mu_d=np.zeros((np.size(MU_D[:,0,0]),180,360))

for i in range(np.size(mu_p[:,0,0])):
    mu_p[i,:,:]=MU_P[i,:,:]*area
    mu_d[i,:,:]=MU_D[i,:,:]*area

#%% Shift dataset
mu_p=np.roll(mu_p,-30)
mu_d=np.roll(mu_d,-30)

#%% Sum over NA (45N - 70N x 270E - 30E)
mu_p1=np.nansum(mu_p[:,135:160,240:],axis=1)
mu_ptot=np.nansum(mu_p1,axis=1)

mu_d1=np.nansum(mu_d[:,135:160,240:],axis=1)
mu_dtot=np.nansum(mu_d1,axis=1)

#%% Moving mean (5 years)
tm=5

mu_pTot = pd.DataFrame(mu_ptot)
rolling_windows = mu_pTot.rolling(tm)
mu_pTOT=np.squeeze(rolling_windows.mean())

mu_dTot = pd.DataFrame(mu_dtot)
rolling_windows = mu_dTot.rolling(tm)
mu_dTOT=np.squeeze(rolling_windows.mean())

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)   
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,mu_pTOT/1e6,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Growth rate (sp) [10$^{-6}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS12_1.png', format='png', dpi=300)

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,mu_dTOT/1e6,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Growth rate (di) [10$^{-6}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS12_2.png', format='png', dpi=300)