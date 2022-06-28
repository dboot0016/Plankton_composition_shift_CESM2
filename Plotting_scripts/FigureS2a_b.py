"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S1a-b

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Amon.co2.gn

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xesmf as xe

#%% CO2
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='co2'    # Variable name

#%% Load in data
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Select and average data
l1=load_var1[var1].mean('lat').mean('lon').mean('plev').compute().squeeze()
l2=load_var2[var1].mean('lat').mean('lon').mean('plev').compute().squeeze()

#%% Concatenate data
VAR=xr.concat([l1,l2],dim='time')

#%% Moving mean (1 year)
tm=12

CO2_0 = pd.DataFrame(VAR)
rolling_windows = CO2_0.rolling(tm)
VAR_gr=np.squeeze(rolling_windows.mean())

#%% Plotting variables
FS=18
t=np.arange(2015+1/24,2101,1/12)
LW=6
a=[2025,2040,2055,2070,2085,2100]

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
im=ax.plot(t,VAR_gr*1e6,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('CO$_2$ [ppm]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS1_1.png', format='png', dpi=300)

#%% Emissions SSP5-8.5 (following ...) 
a=([39152.726,43712.349,55296.583,68775.698,83298.220,100338.606,116805.249,129647.035,130576.239,126287.310])
b=[2015,2020,2030,2040,2050,2060,2070,2080,2090,2100]
c=np.arange(2015,2101,1)

#%% Interpolate data
VAR_gr=np.interp(c,b,a)

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)
LW=6
a=[2025,2040,2055,2070,2085,2100]

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
im=ax.plot(t,VAR_gr/1e3,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Emissions [Pg CO$_2$/yr]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS1_2.png', format='png', dpi=300)