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

suffix = '01'

fn = dir0+'ocean_his_LV4_EPA_winter_000'+suffix+'.nc'

ds = nc.Dataset(fn)

# select wave direction and significant wave height

Dwave = ds['Dwave'][:]
Hwave = ds['Hwave'][:]

nt, ny, nx = Hwave.shape

Dwave_max = np.zeros((nt))
Dwave_min = np.zeros((nt))
Dwave_avg = np.zeros((nt))
Hwave_max = np.zeros((nt))
Hwave_min = np.zeros((nt))
Hwave_avg = np.zeros((nt))

for t in range(nt):
    y_inds,x_inds = np.where(Hwave[t,:,:]!=0)
    Hwave0 = Hwave[t,y_inds,x_inds]
    Dwave0 = Dwave[t,y_inds,x_inds]
    Dwave_max[t] = np.max(Dwave0)
    Dwave_min[t] = np.min(Dwave0)
    Dwave_avg[t] = np.mean(Dwave0)
    Hwave_max[t] = np.max(Hwave0)
    Hwave_min[t] = np.min(Hwave0)
    Hwave_avg[t] = np.mean(Hwave0)

# lat_rho = ds['lat_rho'][:]
# lon_rho = ds['lon_rho'][:]
# mask_rho = ds['mask_rho'][:]
ocean_time = ds['ocean_time'][:]

var_list = 'Dwave_max','Dwave_min','Dwave_avg','Hwave_max','Hwave_min','Hwave_avg','ocean_time'

D = dict()
for var in var_list:
    D[var]=locals()[var]

# outfn = home + 'WQ_data/waves_{:04d}.p'.format(tt)
outfn = home + 'WQ_data/waves_'+suffix+'.p'
pickle.dump(D,open(outfn,'wb'))
