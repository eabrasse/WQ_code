#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
import time
import wqfun
from datetime import datetime, timedelta

home = '/dataSIO/ebrasseale/'


f_list = ['/dataSIO/NADB2017/NADB2017_0/Output/river_tracer_4river_NADB2017_0.nc',\
    '/dataSIO/NADB2018/Input/river_tracer_4river_NADB2018.nc',\
    '/dataSIO/NADB2019/Input/river_tracer_4river_NADB2019_0.nc']

nfiles = len(f_list)

# calculate how many days we'll need
# open first and last files and check time stamps
ds0 = nc.Dataset(f_list[0])
ot0 = ds0['river_time'][:]
dt0 = datetime(1999,1,1) + timedelta(days=ot0[35])

ds_1 = nc.Dataset(f_list[-1])
ot_1 = ds_1['river_time'][:]
dt_1 = datetime(1999,1,1) + timedelta(days=ot_1[-35])

ndays_guess = (dt_1-dt0).days + 1

# build new netCDF file
gf_fn = home+'WQ_data/NADBrivers_daily_gf_NADB2017-2018-2019.nc'

print(f'New netCDF file location: {gf_fn}')
# get rid of the old version, if it exists
try:
    os.remove(gf_fn)
except OSError:
    pass # assume error was because the file did not exist
# ds2 = nc.Dataset(gf_fn, 'w', format='NETCDF3_64BIT_OFFSET')
ds2 = nc.Dataset(gf_fn, 'w', format='NETCDF4')

# Copy dimensions to new netCDF file
for dname, the_dim in ds0.dimensions.items():
    if 'time' in dname:
        # for our new file, time will be ndays long
        ds2.createDimension(dname, ndays_guess)
    else:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

#generate variables
# for all variables that are not time varying, copy them into the new file
tic = time.perf_counter()
var_2gf_list = []
for v_name,varin in ds0.variables.items():
    outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
    
    #EAB 6/7/22: added this line so I don't forget to include metadata next time!
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
    if 'river_time' in varin.dimensions:
        var_2gf_list.append(v_name) # if time-varying, add name to "to godin filter" list
        shape0 = [s for s in varin.shape]
        shape0[0] = ndays_guess #change time index
        outVar_fill = np.zeros((shape0))
        outVar[:] = outVar_fill[:]
    else:
        outVar[:] = varin[:]
toc = time.perf_counter()
init_var_time = f"Initializing variables took {toc-tic:0.4f} seconds"
print(init_var_time)

ds0.close()
ds_1.close()

#different problems because of overlap...

print('Begin filtering')
tic0 = time.perf_counter()
nvars = len(var_2gf_list)
total_days = 0
for f in range(nfiles):
    # if mod(f,5):
    print(f'Working on file {(f+1):} of {nfiles:}'+'\n'+f_list[f])
    ds = nc.Dataset(f_list[f])
    
    rt = ds['river_time'][:]
    
    # if only one file, we'll lose two days, one at the beginning and one at the end
    # if it's the first file of many, we'll lose one day at the beginning
    # if it's the last file of many, we'll lose one day at the end
    # so assume ndays will be at minimum (nt/24 - 2), with up to two days 'restored'
    # how many days 'restored' will be ndays_mod
    ndays_mod = 0
    if f>0: #open previous file UNLESS f=0
        ds0 = nc.Dataset(f_list[f-1])
        rt0 = ds0['river_time'][:]
        r0 = np.argwhere(rt>rt0[-1])[0][0]
        ndays_mod+=1
    else:
        r0 = 0
    if f<nfiles-1: #open following file UNLESS f = nflies-1
        ds1 = nc.Dataset(f_list[f+1])
        rt1 = ds1['river_time'][:]
        r1 = np.argwhere(rt1>rt[-1])[0][0]
        ndays_mod+=1
    else:
        r1=-1
    nt = rt[r0:].shape[0]
    ndays = int(nt/24)-2+ndays_mod

    
    #loop through variables, filtering and saving
    varcount=0
    for var_name in var_2gf_list:
        print(f'   filtering {var_name}...')
        tic = time.perf_counter()
        var = ds[var_name]
        dim = len(var.shape)
        
        # use r0 to trim overlap with previous file
        if dim==1:
            var = var[r0:]
        else:
            var = var[r0:,:]
        dim = len(var.shape)
        
        #godin filtering uses 35 hourly time steps
        # to get one value for the first and last days of each file,
        # we need to include extra hourly time steps from 
        # the files directly before and after.
        # For the first and last files, those won't be available
        if f>0: #don't try to read in preceding file for f=0
            if dim==1: # if one dimensional
                var0 = ds0[var_name][-25:]
            else: # if multidimensional
                var0 = ds0[var_name][-25:,:]
            var = np.insert(var,0,var0,axis=0) # add earlier file data to start of array
        
        if f<nfiles-1: # don't try to read in following file when f = nfiles-1
            if dim==1:
                # use r1 to trim overlap with next file
                var1 = ds1[var_name][r1:(r1+24)]
            else:
                var1 = ds1[var_name][r1:(r1+24),:] 
            var = np.append(var,var1,axis=0) # add later file data to end of array
        
        # pick out daily values
        
        if dim==1: # if one dimensional
            var_gf = wqfun.filt_godin(var)
            ds2[var_name][total_days:total_days+ndays] = var_gf[35:-35:24] #trim leading & trailing nans and subsample 1/day
        else: # if multidimensional
            # filter
            var_gf = wqfun.filt_godin_mat(var)
            ds2[var_name][total_days:total_days+ndays,:] = var_gf[35:-35:24,:] #trim leading & trailing nans and subsample 1/day
        varcount+=1
        toc = time.perf_counter()
        var_gf_time = f"   filtering {var_name} took {toc-tic:0.4f} seconds"
        print(var_gf_time)
    # clean up
    total_days+=ndays
    if f>0:
        ds0.close()
    ds.close()
    if f<nfiles-1:
        ds1.close()

print('Finished!')
toc = time.perf_counter()
total_gf_time = f"Filtering all variables took {toc-tic0:0.4f} seconds"
print(total_gf_time)

dt = [datetime(1999,1,1)+timedelta(days=rt) for rt in ds2['river_time'][:]]
dt_diff = [dt[i]-dt[i-1] for i in range(1,len(dt))]
for diff in dt_diff:
    print(diff)
dt_diff_list = [dt[i] for i in range(0,len(dt)-1) if dt_diff[i]>timedelta(days=1)]
if len(dt_diff_list)==0:
    print('no steps were too long')
else:
    print('some time steps were too long')
    for dtd in dt_diff_list:
        print(dtd)
    
ds2.close()