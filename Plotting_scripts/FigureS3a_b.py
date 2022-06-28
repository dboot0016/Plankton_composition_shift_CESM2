"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S3a-b

Necessary datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Omon.fgco2.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Amon.co2.gn
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Lmon.nbp.gn

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xesmf as xe

#%% Gas exchange
data = '/Users/daan/CESM2_data'   # Location of dataset(s)
var1='fgco2'    # Variable name

#%% Load in and select data
load_var1 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_2.nc')

lv1=load_var1[var1].squeeze()
lv2=load_var2[var1].squeeze()

fgco2_esm=xr.concat([lv1,lv2],dim="time")

#%% Integrate over area in gn grid
load_var1 = xr.open_dataset(f'{data}/area_gn.nc')
area=load_var1['areacello'].squeeze()

fgCO2_esm=(fgco2_esm*area).sum(['nlat','nlon'])

#%% CO2
data = '/Users/daan/CESM2_data'   # Location of dataset(s)
var1='co2'    # Variable name

#%% Load in and select data
load_var1 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_2.nc')

lv1=load_var1[var1].mean('plev').mean('lat').mean('lon').squeeze()
lv2=load_var2[var1].mean('plev').mean('lat').mean('lon').squeeze()

CO2_esm=xr.concat([lv1,lv2],dim="time")

#%% Differentiate w.r.t. time
dCO2_esm=(CO2_esm*44/28.966*5.148e18).differentiate("time")

#%% NBP
data = '/Users/daan/CESM2_data'   # Location of dataset(s)
var1='nbp'    # Variable name

#%% Load in and select data
load_var1 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data}/'+var1+'_esm585_1_2.nc')

lv1=load_var1[var1].squeeze()
lv2=load_var2[var1].squeeze()

nbp_esm=xr.concat([lv1,lv2],dim="time")

#%% Integrating NBP
lon=nbp_esm.lon
lat=nbp_esm.lat

RE = 6.371e6  # [m] Earth radius

dy = 2*np.pi*RE*(lat[1]-lat[0]).values/360                              # grid size in y-direction
dx = 2*np.pi*RE*((lon[1]-lon[0]).values*np.cos(np.deg2rad(lat)))/360    # grid size in x-direction

VAR_lon=(nbp_esm[:,:,:]*dx).sum(['lon'])
NBP_esm=(VAR_lon[:,:]*dy).sum(['lat'])

#%% Determine emissions by summing up
EM_esm=dCO2_esm+NBP_esm+fgCO2_esm

#%% Moving mean (5 years)
tm=60

CO2_0 = pd.DataFrame(dCO2_esm)
rolling_windows = CO2_0.rolling(tm)
dCO2_ESM=np.squeeze(rolling_windows.mean())

CO2_0 = pd.DataFrame(NBP_esm)
rolling_windows = CO2_0.rolling(tm)
NBP_ESM=np.squeeze(rolling_windows.mean())

CO2_0 = pd.DataFrame(fgCO2_esm)
rolling_windows = CO2_0.rolling(tm)
fgCO2_ESM=np.squeeze(rolling_windows.mean())

CO2_0 = pd.DataFrame(EM_esm)
rolling_windows = CO2_0.rolling(tm)
EM_ESM=np.squeeze(rolling_windows.mean())

CO2_0 = pd.DataFrame(CO2_esm)
rolling_windows = CO2_0.rolling(tm)
CO2_ESM=np.squeeze(rolling_windows.mean())

#%% Transform units
EM_ESM=EM_ESM*365*86400*1e-12
fgCO2_ESM=fgCO2_ESM*365*86400*1e-12
dCO2_ESM=dCO2_ESM*365*86400*1e-12
NBP_ESM=NBP_ESM*365*86400*1e-12

#%% Determine relative ratios
r14=np.array(fgCO2_ESM)[:-1]/np.array(EM_ESM)
r13=np.array(NBP_ESM)[:-1]/np.array(EM_ESM)
r12=np.array(dCO2_ESM)[:-1]/np.array(EM_ESM)

R14=np.nancumsum(np.array(fgCO2_ESM)[:-1],axis=0)/np.nancumsum(np.array(EM_ESM),axis=0)
R13=np.nancumsum(np.array(NBP_ESM)[:-1],axis=0)/np.nancumsum(np.array(EM_ESM),axis=0)
R12=np.nancumsum(np.array(dCO2_ESM)[:-1],axis=0)/np.nancumsum(np.array(EM_ESM),axis=0)

#%% Plotting variables
FS=18
t=np.arange(2015+3/24,2101,1/12)
LW=6
a=[2025,2040,2055,2070,2085,2100]

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])

plt.plot(t,np.array(r14)+np.array(r12),linewidth=LW)
plt.plot(t,np.array(r12),linewidth=LW)
plt.plot(t,np.array(r13)+np.array(r14)+np.array(r12),linewidth=LW+3)

ax.fill_between(t, np.array(r12),np.array(r14)+np.array(r12),alpha=0.5)
ax.fill_between(t, 0, np.array(r12),alpha=0.5)
ax.fill_between(t, np.array(r14)+np.array(r12),1,alpha=0.5)
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Ratio of emissions [-]',fontsize=FS)

plt.xlim([2020,2100])
plt.ylim([0,1])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)

#plt.savefig('figureS3_1.png', format='png', dpi=300)

#%% Plotting variables
FS=18
t=np.arange(2015+3/24,2101,1/12)
LW=6
a=[2025,2040,2055,2070,2085,2100]

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])

plt.plot(t,np.array(R14)+np.array(R12),linewidth=LW)
plt.plot(t,np.array(R12),linewidth=LW)
plt.plot(t,np.array(R13)+np.array(R14)+np.array(R12),linewidth=LW)

ax.fill_between(t, np.array(R12),np.array(R14)+np.array(R12),alpha=0.5)
ax.fill_between(t, 0, np.array(R12),alpha=0.5)
ax.fill_between(t, np.array(R14)+np.array(R12),1,alpha=0.5)
plt.grid()
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Cumulative ratio of emissions [-]',fontsize=FS)

plt.xlim([2020,2100])
plt.ylim([0,1])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.legend(['Ocean','Atmosphere','Terrestrial Biosphere'],fontsize=FS-3)

#plt.savefig('figureS3_2.png', format='png', dpi=300)