"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure 3a-f

Necessary processed datasets:
- SCPM_X_Y.txt
   

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xesmf as xe

#%% Load in feedback strength
Datadir='/Users/daan/Documents/PhD/Project 2/SCPM/'

X_Y=np.loadtxt(Datadir+'SCPM_X_Y.txt')
X=-X_Y[:,1]

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)  
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure 4_1
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,np.cumsum(X),linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Feedback strength [ppm]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figure4_1.png', format='png', dpi=300)

#%% Figure 4_2
GtCO2=np.cumsum(X)*1/((1/(44/12))/2.12) #used to convert ppm CO2 to GtCO2 (~7.77)

#2.12: 1 ppm CO2 = 2.12 GtC
#44/12: C --> CO2

fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,GtCO2,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Feedback strength [GtCO$_2$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figure4_2.png', format='png', dpi=300)
