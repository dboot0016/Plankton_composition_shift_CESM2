"""
Effect of plankton composition shifts in the North Atlantic on atmospheric pCO2
Boot, A., von der Heydt, A.S., and Dijkstra, H.A. (2022)

Script for plotting Figure 3a-f

Necessary ESGF datasets:
- CMIP6.C4MIP.NCAR.CESM2.esm-ssp585.r1i1p1f1.Amon.co2.gn

Necessary processed datasets:
- Balk.txt
- Bdic.txt
- SCPM_X_Y.txt

@author: Amber Boot (d.boot@uu.nl)
"""

#%% Load in packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xesmf as xe

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error, r2_score

#%% Load in Biological fluxes
Datadir='/Users/daan/Documents/PhD/Project 2/SCPM/'

Balk=np.loadtxt(Datadir+'Balk.txt')
Bdic=np.loadtxt(Datadir+'Bdic.txt')

BALK=Balk[4,:]
BDIC=Bdic[4,:]

#%% Load in CO2 data
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='co2'    # Variable name

#%% Load in dataset
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1_2.nc')

#%% Average dataset + select data
l1=load_var1[var1].mean('lat').mean('lon').mean('plev').compute().squeeze()
l2=load_var2[var1].mean('lat').mean('lon').mean('plev').compute().squeeze()

#%% Concatenate data
VAR=xr.concat([l1,l2],dim='time')

#%% Take annual averages of CO2
co2=np.zeros((86,))

for i in range(86):
    co2[i]=np.mean(VAR[i*12:(i+1)*12])

#%% Define a function to be fitted
CO2=co2*1e6
def func(CO2, a, b):
    return a*np.log(CO2)+b

#%% Fit alkalinity flux
from scipy.optimize import curve_fit
popt, pcov = curve_fit(func, CO2, BALK)

#%% Fit DIC flux
from scipy.optimize import curve_fit
popt2, pcov2 = curve_fit(func, CO2, BDIC)

#%% Determine R2
residuals_A = BALK - func(CO2, *popt)
residuals_C = BDIC - func(CO2, *popt2)

ss_res_A=np.sum(residuals_A**2)
ss_res_C=np.sum(residuals_C**2)

ss_tot_A=np.sum((BALK-np.mean(BALK))**2)
ss_tot_C=np.sum((BDIC-np.mean(BDIC))**2)

R2_A=1-(ss_res_A/ss_tot_A)
R2_C=1-(ss_res_C/ss_tot_C)

#%% Create fit dataset
CO2_fit=np.linspace(350,1200,150)
balk=popt[0]*np.log(CO2_fit)+popt[1]
bdic=popt2[0]*np.log(CO2_fit)+popt2[1]

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)  
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure 3_1
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,BALK,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Alkalinity flux [mol m$^{-3}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figure3_1.png', format='png', dpi=300)

#%% Figure 3_2
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,BDIC,linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('DIC flux [mol C m$^{-3}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
#plt.savefig('figure3_2.png', format='png', dpi=300)

#%% Plotting variables
FS=18
a=[400,550,700,850,1000]
LW=6
MS=40 

#%% Figure 3_3
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.scatter(CO2,BALK,s=MS)
ax.plot(CO2_fit,balk,color='red',linewidth=LW)
plt.xlabel('CO$_2$ [ppm]',fontsize=FS)
plt.ylabel('Alkalinity flux [mol m$^{-3}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([350,1100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.legend(['Fit (R$^2$ = '+np.array2string(R2_A, formatter={'float_kind':lambda x: "%.2f" % x})+')','CESM data'],fontsize=FS-3)
#plt.savefig('figure3_3.png', format='png', dpi=300)

#%% Figure 3_4
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.scatter(CO2,BDIC,s=MS)
ax.plot(CO2_fit,bdic,color='red',linewidth=LW)
plt.xlabel('CO$_2$ [ppm]',fontsize=FS)
plt.ylabel('DIC flux [mol C m$^{-3}$ s$^{-1}$]',fontsize=FS)
plt.grid()
plt.xlim([350,1100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.legend(['Fit (R$^2$ = '+np.array2string(R2_C, formatter={'float_kind':lambda x: "%.2f" % x})+')','CESM data'],fontsize=FS-3)
#plt.savefig('figure3_4.png', format='png', dpi=300)

#%% Load in model mismatch + feedback strength
Datadir='/Users/daan/Documents/PhD/Project 2/SCPM/'

X_Y=np.loadtxt(Datadir+'SCPM_X_Y.txt')
X=X_Y[:,1]
Y=X_Y[:,0]

#%% Plotting variables
FS=18
t=np.arange(2015,2101,1)  
a=[2025,2040,2055,2070,2085,2100]
LW=6

#%% Figure 3_5
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,Y,linewidth=LW)
ax.plot(t,-X,color='red',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('$\Delta$CO$_2$/$\Delta$t [ppm/yr]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.ylim([-0.75,4.5])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.legend(['Uncaptured dynamics','Feedback strength'],fontsize=FS-3)
#plt.savefig('figure3_5.png', format='png', dpi=300)

#%% Figure 3_6
fig = plt.figure(figsize=(7, 4.5))

ax = fig.add_axes([0.25,0.15,0.67,0.8])
ax.plot(t,np.cumsum(Y)/(co2*1e6),linewidth=LW)
ax.plot(t,-np.cumsum(X)/(co2*1e6),color='red',linewidth=LW)
plt.xlabel('Time',fontsize=FS)
plt.ylabel('Relative cumulative $\Delta$CO$_2$ [-]',fontsize=FS)
plt.grid()
plt.xlim([2020,2100])
plt.xticks(a,fontsize=FS-3)
plt.yticks(fontsize=FS-3)
plt.tight_layout()
plt.legend(['Uncaptured dynamics','Feedback strength'],fontsize=FS-3)
#plt.savefig('figure3_6.png', format='png', dpi=300)