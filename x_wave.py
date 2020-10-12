#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import pickle
import numpy as np
import netCDF4 as nc

home = '/home/ebrasseale/'

dir0 = '/home/x1wu/SDTJRE_EPA/LV4_RUNFILES/Run2017/Couple/'

fn = dir0+'ocean_his_LV4_EPA_winter_00001.nc'

ds = nc.Dataset(fn)

# select wave direction and significant wave height

# at first time step 
tt = 0

Dwave = ds['Dwave'][tt,:]
Hwave = ds['Hwave'][tt,:]
lat_rho = ds['lat_rho'][:]
lon_rho = ds['lon_rho'][:]
mask_rho = ds['mask_rho'][:]
ocean_time = ds['ocean_time'][tt]

var_list = 'Dwave','Hwave','lat_rho','lon_rho','mask_rho','ocean_time'

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/waves_{:04d}.p'.format(tt)
pickle.dump(D,open(outfn,'wb'))
