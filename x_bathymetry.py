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

home = '/data0/ebrasseale/'

dir0 = '/data0/NADB2017/NADB2017_0_NEW/'

fn = dir0+'ocean_his_NADB_0_new_00001.nc'

ds = nc.Dataset(fn)

# lon_rho = ds['lon_rho'][:]
# lat_rho = ds['lat_rho'][:]
# mask_rho = ds['mask_rho'][:]
# h = ds['h'][:]
D={}
var_list = ['lon_rho','lat_rho','mask_rho','h']
for varname in var_list:
   D[varname] = ds[varname][:] 

outfn = home + 'WQ_data/NADB2017_bathymetry.p'
pickle.dump(D,open(outfn,'wb'))

