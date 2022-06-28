"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure S14a-f

Necessary processed datasets:
- AMOC_esm585_1.nc

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np

#%% AMOC
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='__xarray_dataarray_variable__'    # Variable name

#%% load in data + select
load_var1 = xr.open_dataset(f'{data1}/AMOC_esm585_1.nc')
l1=load_var1[var1].compute().squeeze()

#%% Determine maximum
VAR=np.zeros((1032,))

for i in range(1032):
    VAR[i]=max(l1[i,:])

#%% Turn into xarray
VAR=VAR*xr.ones_like(l1[:,0].squeeze())

#%% Take 5 year running mean
VAR_gr=VAR.rolling(time=60, center=False).mean()

#%% Plotting variables
FS=18
t=np.arange(2015+1/24,2101,1/12)  
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,VAR_gr/1e6,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('AMOC (26.5$^{\circ}$N) [Sv]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figureS7.png', format='png', dpi=300)
