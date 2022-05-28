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

dir0 = '/dataSIO/NADB2017/NADB2017_0_NEW/'
dir1 = '/dataSIO/NADB2018/'
dir2 = '/dataSIO/NADB2019/'

f_list = []
for my_dir in [dir0,dir1,dir2]:
    f_list0 = os.listdir(my_dir)
    f_list0.sort()
    f_list0 = [my_dir+x for x in f_list0 if x[:9]=='ocean_his']
    f_list.extend(f_list0)

testing=False
if testing:
    f_list = f_list[:4]
    print('TESTING!!!')

nfiles = len(f_list)

# calculate how many days we'll need
# open first and last files and check time stamps
ds0 = nc.Dataset(f_list[0])
ot0 = ds0['ocean_time'][:]
dt0 = datetime(1999,1,1) + timedelta(seconds=ot0[35])

ds1 = nc.Dataset(f_list[-1])
ot1 = ds1['ocean_time'][:]
dt1 = datetime(1999,1,1) + timedelta(seconds=ot1[-35])

ndays_guess = (dt1-dt0).days + 1

# build new netCDF file
gf_fn = home+'WQ_data/ocean_daily_gf_NADB2017-2018-2019.nc'

print(f'New netCDF file location: {gf_fn}')
# get rid of the old version, if it exists
try:
    os.remove(gf_fn)
except OSError:
    pass # assume error was because the file did not exist
# ds2 = nc.Dataset(gf_fn, 'w', format='NETCDF3_64BIT_OFFSET')
ds2 = nc.Dataset(gf_fn, 'w', format='NETCDF4_CLASSIC')

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
    if 'ocean_time' in varin.dimensions:
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

#close to clean up, want to make sure we know exactly what's open
ds0.close()
ds1.close()

print('Begin filtering')
tic0 = time.perf_counter()
nvars = len(var_2gf_list)
total_days = 0
for f in range(nfiles):
    # if mod(f,5):
    print(f'Working on file {(f+1):} of {nfiles:}'+'\n'+f_list[f])
    ds = nc.Dataset(f_list[f])
    
    # if only one file, we'll lose two days, one at the beginning and one at the end
    # if it's the first file of many, we'll lose one day at the beginning
    # if it's the last file of many, we'll lose one day at the end
    # so assume ndays will be at minimum (nt/24 - 2), with up to two days 'restored'
    # how many days 'restored' will be ndays_mod
    ndays_mod = 0
    if f>0: #open previous file UNLESS f=0
        ds0 = nc.Dataset(f_list[f-1])
        ndays_mod+=1
    if f<nfiles-1: #open following file UNLESS f = nflies-1
        ds1 = nc.Dataset(f_list[f+1])
        ndays_mod+=1
    nt = ds['ocean_time'][:].shape[0]
    ndays = int(nt/24)-2+ndays_mod

    
    #loop through variables, filtering and saving
    varcount=0
    for var_name in var_2gf_list:
        print(f'   filtering {var_name}...')
        tic = time.perf_counter()
        var = ds[var_name]
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
                var1 = ds1[var_name][:24]
            else:
                var1 = ds1[var_name][:24,:] 
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
ds2.close()
print()