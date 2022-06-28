"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S14a-f

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.uo.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.vo.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Oyr.phypico.gn

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xesmf as xe

#%% u-velocity
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='uo'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select data (top 150m)
l1=load_var1[var1][:,:15,:,:].compute().squeeze()
l2=load_var2[var1][:,:15,:,:].compute().squeeze()

#%% Concatenate data
VAR=(xr.concat([l1,l2],dim='time'))

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
uo=np.roll(VAR_gr,-210)

#%% v-velocity
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='vo'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select data (top 150m)
l1=load_var1[var1][:,:15,:,:].compute().squeeze()
l2=load_var2[var1][:,:15,:,:].compute().squeeze()

#%% Concatenate data
VAR=(xr.concat([l1,l2],dim='time'))

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
vo=np.roll(VAR_gr,-210)

#%% Empty variables
l1=[]
l2=[]
VAR=[]
VAR_gr=[]

#%% Select NA boundaries (45N - 70N x 270E - 30E)
uo_W=uo[:,:,135:160,240]
uo_E=uo[:,:,135:160,-30]

vo_N=vo[:,:,160,240:-30]
vo_S=vo[:,:,135,240:-30]

#%% Small phytoplankton biomass
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='phypico'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select data
l1=load_var1[var1].compute().squeeze()
l2=load_var2[var1].compute().squeeze()

#%% Concatenate data
VAR=(xr.concat([l1,l2],dim='time'))

#%% Regrid data
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(VAR, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(VAR)
sp=np.roll(VAR_gr,-210)

#%% Select NA boundaries (45N - 70N x 270E - 30E)
sp_W=sp[:,:,135:160,240]
sp_E=sp[:,:,135:160,-30]

sp_N=sp[:,:,160,240:-30]
sp_S=sp[:,:,135,240:-30]

#%% Annual means velocity data
uO_E=np.zeros((86,15,25))
uO_W=np.zeros((86,15,25))

vO_N=np.zeros((86,15,90))
vO_S=np.zeros((86,15,90))

for i in range(86):
    uO_E[i,:,:]=np.mean(uo_E[i*12:(i+1)*12,:,:],axis=0)
    uO_W[i,:,:]=np.mean(uo_W[i*12:(i+1)*12,:,:],axis=0)
    vO_N[i,:,:]=np.mean(vo_N[i*12:(i+1)*12,:,:],axis=0)
    vO_S[i,:,:]=np.mean(vo_S[i*12:(i+1)*12,:,:],axis=0)
    
#%% Determine fluxes (x10 because layer thickness is 10m)
flux_E=np.nansum(sp_E*uO_E,axis=1)*10
flux_W=np.nansum(sp_W*uO_W,axis=1)*10
flux_N=np.nansum(sp_N*vO_N,axis=1)*10
flux_S=np.nansum(sp_S*vO_S,axis=1)*10

#%%
flux_E1=np.nansum(np.mean(sp_E,axis=0)*uO_E,axis=1)*10
flux_W1=np.nansum(np.mean(sp_W,axis=0)*uO_W,axis=1)*10
flux_N1=np.nansum(np.mean(sp_N,axis=0)*vO_N,axis=1)*10
flux_S1=np.nansum(np.mean(sp_S,axis=0)*vO_S,axis=1)*10

flux_E2=np.nansum(sp_E*np.mean(uO_E,axis=0),axis=1)*10
flux_W2=np.nansum(sp_W*np.mean(uO_W,axis=0),axis=1)*10
flux_N2=np.nansum(sp_N*np.mean(vO_N,axis=0),axis=1)*10
flux_S2=np.nansum(sp_S*np.mean(vO_S,axis=0),axis=1)*10

#%%
U_E=np.nansum(uO_E,axis=1)*10
U_W=np.nansum(uO_W,axis=1)*10
V_N=np.nansum(vO_N,axis=1)*10
V_S=np.nansum(vO_S,axis=1)*10

#%%
RE = 6.371e6  # [m] Earth radius

dy = 2*np.pi*RE*(1)/360                              # grid size in y-direction
dx1 = 2*np.pi*RE*(1*np.cos(np.deg2rad(70)))/360    # grid size in x-direction
dx2 = 2*np.pi*RE*(1*np.cos(np.deg2rad(45)))/360    # grid size in x-direction

F_E=np.nansum(flux_E,axis=1)*dy
F_W=np.nansum(flux_W,axis=1)*dy
F_N=np.nansum(flux_N,axis=1)*dx1
F_S=np.nansum(flux_S,axis=1)*dx2

#%%
F_E1=np.nansum(flux_E1,axis=1)*dy
F_W1=np.nansum(flux_W1,axis=1)*dy
F_N1=np.nansum(flux_N1,axis=1)*dx1
F_S1=np.nansum(flux_S1,axis=1)*dx2

F_E2=np.nansum(flux_E2,axis=1)*dy
F_W2=np.nansum(flux_W2,axis=1)*dy
F_N2=np.nansum(flux_N2,axis=1)*dx1
F_S2=np.nansum(flux_S2,axis=1)*dx2

#%% Running mean (5 years) + balance
tm=5

Var_tot = pd.DataFrame(F_E)
rolling_windows = Var_tot.rolling(tm)
fe=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_W)
rolling_windows = Var_tot.rolling(tm)
fw=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_N)
rolling_windows = Var_tot.rolling(tm)
fn=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_S)
rolling_windows = Var_tot.rolling(tm)
fs=np.squeeze(rolling_windows.mean())

ft=-fe+fw-fn+fs

#%% Running mean (5 years) + balance
Var_tot = pd.DataFrame(F_E1)
rolling_windows = Var_tot.rolling(tm)
fe1=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_W1)
rolling_windows = Var_tot.rolling(tm)
fw1=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_N1)
rolling_windows = Var_tot.rolling(tm)
fn1=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_S1)
rolling_windows = Var_tot.rolling(tm)
fs1=np.squeeze(rolling_windows.mean())

ft1=-fe1+fw1-fn1+fs1

Var_tot = pd.DataFrame(F_E2)
rolling_windows = Var_tot.rolling(tm)
fe2=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_W2)
rolling_windows = Var_tot.rolling(tm)
fw2=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_N2)
rolling_windows = Var_tot.rolling(tm)
fn2=np.squeeze(rolling_windows.mean())

Var_tot = pd.DataFrame(F_S2)
rolling_windows = Var_tot.rolling(tm)
fs2=np.squeeze(rolling_windows.mean())

ft2=-fe2+fw2-fn2+fs2

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)  
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure 1
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,-fe*365*86400/1e10,linewidth=LW)
ax.plot(t,-fe1*365*86400/1e10,color='red',linewidth=LW)
ax.plot(t,-fe2*365*86400/1e10,color='black',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Flux east [10$^{10}$ mol C yr$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.ylim([-1.65,-0.55])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.legend(['$\overrightarrow{V}$ P$_{sp}$','$\overrightarrow{V}$ $\langle$P$_{sp}$⟩','$\langle\overrightarrow{V}$⟩ P$_{sp}$'],fontsize=FS-3,ncol=3)
plt.savefig('figureS14_1.png', format='png', dpi=300)

#%% Figure 2
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,fw*365*86400/1e10,linewidth=LW)
ax.plot(t,fw1*365*86400/1e10,color='red',linewidth=LW)
ax.plot(t,fw2*365*86400/1e10,color='black',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Flux west [10$^{10}$ mol C yr$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.ylim([0,0.21])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.legend(['$\overrightarrow{V}$ P$_{sp}$','$\overrightarrow{V}$ $\langle$P$_{sp}$⟩','$\langle\overrightarrow{V}$⟩ P$_{sp}$'],fontsize=FS-3,ncol=3)
plt.savefig('figureS14_2.png', format='png', dpi=300)

#%% Figure 3
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,-fn*365*86400/1e10,linewidth=LW)
ax.plot(t,-fn1*365*86400/1e10,color='red',linewidth=LW)
ax.plot(t,-fn2*365*86400/1e10,color='black',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Flux north [10$^{10}$ mol C yr$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.savefig('figureS14_3.png', format='png', dpi=300)

#%% Figure 4
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,fs*365*86400/1e10,linewidth=LW)
ax.plot(t,fs1*365*86400/1e10,color='red',linewidth=LW)
ax.plot(t,fs2*365*86400/1e10,color='black',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Flux south [10$^{10}$ mol C yr$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.savefig('figureS14_4.png', format='png', dpi=300)

#%% Figure 5
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,ft*365*86400/1e10,linewidth=LW)
ax.plot(t,ft1*365*86400/1e10,color='red',linewidth=LW)
ax.plot(t,ft2*365*86400/1e10,color='black',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Balance [10$^{10}$ mol C yr$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.savefig('figureS14_5.png', format='png', dpi=300)

#%% Integrate over depth and select NA (45N - 70N x 270E - 30E)
SP=np.nansum(sp[:,:,135:160,240:-30],axis=1)*10

#%% Load in area data
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()
Area=area[135:160,240:-30]

#%% Multiply with area + sum over region
VAR_sp=np.zeros((np.shape(SP)))

for i in range(np.size(SP[:,0,0])):
    VAR_sp[i,:,:]=SP[i,:,:]*Area
    
VAR_SP=np.nansum(np.nansum(VAR_sp,axis=1),axis=1)

#%%
Var_tot = pd.DataFrame(VAR_SP)
rolling_windows = Var_tot.rolling(tm)
var=np.squeeze(rolling_windows.mean())

#%% Figure 6
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,var/1e10,linewidth=LW)

plt.xlabel('Time',fontsize=FS)
plt.ylabel('Biomass (sp) [10$^{10}$ mol C]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.savefig('figureS14_6sp.png', format='png', dpi=300)