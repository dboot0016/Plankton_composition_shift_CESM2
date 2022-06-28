"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S15a-f

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limndiat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limfediat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.limirrdiat.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.phydiat.gn    
    
Necessary processed datasets:
- limsidiat_esm585_1.nc
- limpdiat_esm585_1.nc

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
var1='limndiat'                  # Nitrogen limitation diatoms
var2='limfediat'                 # Iron limitation diatoms
var3='limirrdiat'                # Light limitation diatoms
var4='limsidiat'                 # Silicate limitation diatoms
var5='limpdiat'                  # Nitrogen limitation small phytoplankton
            
#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')
load_var3 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_1.nc')
load_var4 = xr.open_dataset(f'{data1}/'+var2+'_esm585_1_2.nc')
load_var5 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_1.nc')
load_var6 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_2.nc')
load_var7 = xr.open_dataset(f'{data1}/'+var4+'_esm585_1.nc')
load_var8 = xr.open_dataset(f'{data1}/'+var5+'_esm585_1.nc')

#%%
l1=load_var1[var1].compute().squeeze()
l2=load_var2[var1].compute().squeeze()
l3=load_var3[var2].compute().squeeze()
l4=load_var4[var2].compute().squeeze()
l5=load_var5[var3].compute().squeeze()
l6=load_var6[var3].compute().squeeze()

#%%
varx='__xarray_dataarray_variable__'
l7=load_var7[varx].compute().squeeze()
l8=load_var8[varx].compute().squeeze()

#%%
limN_diat=xr.concat([l1,l2],dim='time')
limF_diat=xr.concat([l3,l4],dim='time')
limI_diat=xr.concat([l5,l6],dim='time')
limS_diat=l7
limP_diat=l8

#%%
limf_Di=np.array(limF_diat)
limn_Di=np.array(limN_diat)
lims_di=np.array(limS_diat)
limp_di=np.array(limP_diat)
limI_Diat=np.array(limI_diat)

#%%
limf_di=np.zeros((86,384,320))
limn_di=np.zeros((86,384,320))
limi_diat=np.zeros((86,384,320))

for i in range(86):
    limf_di[i,:,:]=np.mean(limf_Di[i*12:(i+1)*12,:,:],axis=0)
    limn_di[i,:,:]=np.mean(limn_Di[i*12:(i+1)*12,:,:],axis=0)
    limi_diat[i,:,:]=np.mean(limI_Diat[i*12:(i+1)*12,:,:],axis=0)
    
#%%
limdi=np.zeros((86,384,320))

for i in range(384):
    for j in range(320):
        for k in range(86):
            limdi[k,i,j]=min([limn_di[k,i,j],limp_di[k,i,j],limf_di[k,i,j],lims_di[k,i,j]])
               
#%%
for i in range (384):
    for j in range(320):
        if np.isnan(limn_di[0,i,j]):
            limdi[:,i,j]='nan'

#%%
co_lim=limdi*limi_diat    
co_lim=co_lim*xr.ones_like(l7)

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(co_lim, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(co_lim)
var_gr=np.roll(VAR_gr,-180)

#%% Load in area data
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%% Multiply with area
VAR_gr=np.zeros((np.size(var_gr[:,0,0]),180,360))

for i in range(np.size(var_gr[:,0,0])):
    VAR_gr[i,:,:]=var_gr[i,:,:]*area

#%% Determine where nans are
nans=var_gr[0,:,:]/var_gr[0,:,:]

#%% Shift datasets
VAR_gr=np.roll(VAR_gr,-30)
area_gr=np.roll(area*nans,-30)

#%% Integrate over NA (45N - 70N x 270E - 30E)
VAR_1=np.nansum(VAR_gr[:,135:160,240:],axis=1)
VAR_tot=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[135:160,240:],axis=1)
area_tot=np.nansum(area_1,axis=0)

#%% Running mean + average over region
tm=5

Var_tot = pd.DataFrame(VAR_tot)
rolling_windows = Var_tot.rolling(tm)
var_tot=np.squeeze(rolling_windows.mean())

Co_lim=var_tot/area_tot

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)
a=[2025,2040,2055,2070,2085,2100]
b=[0.170,0.175,0.180,0.185]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,Co_lim,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Co-imitation (di) [-]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(b,fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS10_4.png', format='png', dpi=300)

#%% Light limitation diatoms
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='limirrdiat'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select data
l1=load_var1[var1].compute().squeeze()
l2=load_var2[var1].compute().squeeze()

#%% Concatenate in time
VAR=(xr.concat([l1,l2],dim='time'))

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
var_gr=np.roll(VAR_gr,-180)

#%% Load in area data
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%% Multiply with area
VAR_gr=np.zeros((np.size(var_gr[:,0,0]),180,360))

for i in range(np.size(var_gr[:,0,0])):
    VAR_gr[i,:,:]=var_gr[i,:,:]*area

#%% Determine where nans are
nans=var_gr[0,:,:]/var_gr[0,:,:]

#%% Shift datasets
VAR_gr=np.roll(VAR_gr,-30)
area_gr=np.roll(area*nans,-30)

#%% Integrate over NA (45N - 70N x 270E - 30E)
VAR_1=np.nansum(VAR_gr[:,135:160,240:],axis=1)
VAR_tot=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[135:160,240:],axis=1)
area_tot=np.nansum(area_1,axis=0)

#%% Running mean + average over region
tm=60

Var_tot = pd.DataFrame(VAR_tot)
rolling_windows = Var_tot.rolling(tm)
var_tot=np.squeeze(rolling_windows.mean())

varL_tot=var_tot/area_tot

#%% Plotting variables
FS=18
t=np.arange(2015+1/24,2101,1/12)
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,varL_tot,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Light limitation (di) [-]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()

#plt.savefig('figureS10_1.png', format='png', dpi=300)

#%% Nitrogen limitation diatoms
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='limndiat'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select data
l1=load_var1[var1].compute().squeeze()
l2=load_var2[var1].compute().squeeze()

#%% Concatenate
VAR=(xr.concat([l1,l2],dim='time'))

#%% Regrid
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
var_gr=np.roll(VAR_gr,-180)

#%% Load in area data
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%% Multiply with area
VAR_gr=np.zeros((np.size(var_gr[:,0,0]),180,360))

for i in range(np.size(var_gr[:,0,0])):
    VAR_gr[i,:,:]=var_gr[i,:,:]*area

#%% Determine where nans
nans=var_gr[0,:,:]/var_gr[0,:,:]

#%% Shift dataset 
VAR_gr=np.roll(VAR_gr,-30)
area_gr=np.roll(area*nans,-30)

#%% Integrate over NA (45N - 70N x 270E - 30E)
VAR_1=np.nansum(VAR_gr[:,135:160,240:],axis=1)
VAR_tot=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[135:160,240:],axis=1)
area_tot=np.nansum(area_1,axis=0)

#%% Running mean and average over area
tm=60

Var_tot = pd.DataFrame(VAR_tot)
rolling_windows = Var_tot.rolling(tm)
var_tot=np.squeeze(rolling_windows.mean())

varN_tot=var_tot/area_tot

#%% Plotting variables
FS=18
t=np.arange(2015+1/24,2101,1/12)   
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,varN_tot,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('N limitation (di) [-]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS10_2.png', format='png', dpi=300)

#%% Silicate limitation diatoms
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='limsidiat'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1.nc')
varx='__xarray_dataarray_variable__'
VAR=load_var1[varx].compute().squeeze()

#%% Regrid
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
var_gr=np.roll(VAR_gr,-180)

#%% Load in area data
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

#%% Multiply with area
VAR_gr=np.zeros((np.size(var_gr[:,0,0]),180,360))

for i in range(np.size(var_gr[:,0,0])):
    VAR_gr[i,:,:]=var_gr[i,:,:]*area

#%% Shift dataset 
VAR_gr=np.roll(VAR_gr,-30)
area_gr=np.roll(area*nans,-30)

#%% Integrate over NA (45N - 70N x 270E - 30E)
VAR_1=np.nansum(VAR_gr[:,135:160,240:],axis=1)
VAR_tot=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[135:160,240:],axis=1)
area_tot=np.nansum(area_1,axis=0)

#%% Running mean and average over area
tm=5

Var_tot = pd.DataFrame(VAR_tot)
rolling_windows = Var_tot.rolling(tm)
var_tot=np.squeeze(rolling_windows.mean())

varS_tot=var_tot/area_tot

#%% Plotting variables
FS=18
t2=np.arange(2015,2101,1)   
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t2,varS_tot,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Si limitation (di) [-]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS10_3.png', format='png', dpi=300)