import xarray as xr 
import numpy as np
import xesmf as xe

#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='bddtalk'    # Variable name

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_2.nc')

#%%
lb=load_var1[var1][:,:10,:,:].mean('lev').compute().squeeze()
le=load_var2[var1][:,:10,:,:].mean('lev').compute().squeeze()

#%%
l2=xr.concat([lb,le],dim='time')

#%%
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(l2, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(l2)
L2=np.roll(VAR_gr,0)

#%%
Lat1=130
Lat3=150

Lon1=90
Lon3=210

#%%
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

area=np.roll(area,-180)

#%%
LL2=area*L2

#%%
nans=L2[0,:,:]/L2[0,:,:]
area_gr=area*nans

#%%
VAR_1=np.nansum(LL2[:,Lat1:Lat3,Lon1:Lon3],axis=1)
Val5=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[Lat1:Lat3,Lon1:Lon3],axis=1)
A2=np.nansum(area_1,axis=0)

#%%
BAlk=Val5/A2*86400*365

#%%
#np.savetxt('Balk.txt',BAlk)
#%%
data1='/Users/daan/CESM2_data'   # Location of dataset(s)
var1='bddtdic'    # Variable name

#%%
load_var1 = xr.open_dataset(f'{data1}/'+var1+'_esm585_1.nc')
load_var2 = xr.open_dataset(f'{data1}/'+var1+'_esm585_2.nc')

#%%
lb=load_var1[var1][:,:10,:,:].mean('lev').compute().squeeze()
le=load_var2[var1][:,:10,:,:].mean('lev').compute().squeeze()

#%%
l2=xr.concat([lb,le],dim='time')

#%%
ds_out = xe.util.grid_global(1, 1)
regridder = xe.Regridder(l2, ds_out, 'bilinear',periodic=True)
VAR_gr = regridder(l2)
L2=np.roll(VAR_gr,0)

#%%
Lat1=130
Lat3=150

Lon1=90
Lon3=210

#%%
load_var1 = xr.open_dataset(f'{data1}/area_gr.nc')
lat=load_var1['lat'].load()
lon=load_var1['lon'].load()
area=load_var1['areacello'].load()

area=np.roll(area,-180)

#%%
LL2=area*L2

#%%
nans=L2[0,:,:]/L2[0,:,:]
area_gr=area*nans

#%%
VAR_1=np.nansum(LL2[:,Lat1:Lat3,Lon1:Lon3],axis=1)
Val5=np.nansum(VAR_1,axis=1)

area_1=np.nansum(area_gr[Lat1:Lat3,Lon1:Lon3],axis=1)
A2=np.nansum(area_1,axis=0)

#%%
BDIC=Val5/A2*86400*365

#%%
#np.savetxt('BDIC.txt',BDIC)