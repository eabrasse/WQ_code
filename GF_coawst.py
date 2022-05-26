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

testing=True
if testing:
    f_list = f_list[:2]
    print('TESTING!!!')

nfiles = len(f_list)
ndays_guess = nfiles*10 # not every file has ten days, but most do

# build new netCDF file
gf_fn = home+'WQ_data/ocean_daily_gf_NADB2017-2018-2019.nc'

print(f'New netCDF file location: {gf_fn}')
# get rid of the old version, if it exists
try:
    os.remove(gf_fn)
except OSError:
    pass # assume error was because the file did not exist
ds2 = nc.Dataset(gf_fn, 'w', format='NETCDF3_64BIT_OFFSET')

#open first file
ds = nc.Dataset(f_list[0])

# Copy dimensions to new netCDF file
for dname, the_dim in ds.dimensions.items():
    if 'time' in dname:
        # for our new file, time will be ndays long
        ds2.createDimension(dname, ndays_guess)
    else:
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

#generate variables
# for all variables that are not time varying, copy them into the new file
var_2gf_list = []
for v_name,varin in ds.variables.items():
    outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
    if 'ocean_time' in varin.dimensions:
        var_2gf_list.append(v_name) # if time-varying, add name to "to godin filter" list
        shape = [s for s in varin.shape]
        shape[0] = ndays_guess #change time index

        outVar[:] = np.zeros((shape))
    else:
        outVar[:] = varin[:]

nvars = len(var_2gf_list)
total_days = 0
for f in range(nfiles):
    # if mod(f,5):
    print(f'Working on file {(f+1):} of {nfiles:}'+'\n'+f_list[f])
    ds = nc.Dataset(f_list[f])
    if f>0: #open previous file UNLESS f=0
        ds0 = nc.Dataset(f_list[f-1])
    if f<nfiles-1: #open following file UNLESS f = nflies-1
        ds1 = nc.Dataset(f_list[f+1])
    nt = ds['ocean_time'][:].shape[0]
    ndays = int(nt/24)
    
    #loop through variables, filtering and saving
    varcount=0
    for var_name in var_2gf_list:
        print(f'   filtering {var_name}...')
        var = ds[var_name]
        dim = len(var.shape)
        
        #godin filtering uses 35 hourly time steps
        # to get one value for the first and last days of each file,
        # we need to include extra hourly time steps from 
        # the files directly before and after.
        # For the first and last files, those won't be available
        if f>0: #don't try to read in preceding file for f=0
            if dim==1: # if one dimensional
                var0 = ds0[var_name][-24:]
            else: # if multidimensional
                var0 = ds0[var_name][-24:,:]
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
    # clean up
    total_days+=ndays
    if f>0:
        ds0.close()
        
print('Finished!')
ds.close() # since no new ds1 is loaded for the last day, ds1 should refer to the same file as ds
ds2.close()
print()