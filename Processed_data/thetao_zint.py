
import xarray as xr 
import numpy as np

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var3='thetao'    # Variable name

#%%
load_var5 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_1.nc')
load_var6 = xr.open_dataset(f'{data1}/'+var3+'_esm585_1_2.nc')

#%%
l5=load_var5[var3][:,:10,:,:].sum('lev').compute().squeeze()/10
l6=load_var6[var3][:,:10,:,:].sum('lev').compute().squeeze()/10

L3=xr.concat([l5,l6],dim='time')

#%%
L5=np.zeros((86,384,320))

for i in range(86):
    L5[i,:,:]=L3[i*12:(i+1)*12,:,:].mean('time')

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var4='limsidiat'    # Variable name

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var4+'_esm585_1.nc')
L3=load_var1[var4].squeeze()

#%%
L5=L5*xr.ones_like(L3)
                           
#%%
#L5.to_netcdf('thetao_esm585_zint_1.nc')
