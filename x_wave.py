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

# Dwave_max = np.zeros((nt))
# Dwave_min = np.zeros((nt))
Dwave_avg = np.zeros((nt))
Dwave_10 = np.zeros((nt))
Dwave_25 = np.zeros((nt))
Dwave_75 = np.zeros((nt))
Dwave_90 = np.zeros((nt))
# Dwave_std = np.zeros((nt))
# Hwave_max = np.zeros((nt))
# Hwave_min = np.zeros((nt))
Hwave_avg = np.zeros((nt))
Hwave_10 = np.zeros((nt))
Hwave_25 = np.zeros((nt))
Hwave_75 = np.zeros((nt))
Hwave_90 = np.zeros((nt))
# Hwave_std = np.zeros((nt))

for t in range(nt):
    y_inds,x_inds = np.where(Hwave[t,:,:]!=0)
    Hwave0 = Hwave[t,y_inds,x_inds]
    Hwave0 = np.ma.filled(Hwave0,np.nan)
    Dwave0 = Dwave[t,y_inds,x_inds]
    Dwave0 = np.ma.filled(Dwave0,np.nan)
    
    Dwave_avg[t] = np.nanmean(Dwave0)
    Dwave_10[t] = np.nanpercentile(Dwave0,10)
    Dwave_25[t] = np.nanpercentile(Dwave0,25)
    Dwave_75[t] = np.nanpercentile(Dwave0,75)
    Dwave_90[t] = np.nanpercentile(Dwave0,90)
    Hwave_avg[t] = np.nanmean(Hwave0)
    Hwave_10[t] = np.nanpercentile(Hwave0,10)
    Hwave_25[t] = np.nanpercentile(Hwave0,25)
    Hwave_75[t] = np.nanpercentile(Hwave0,75)
    Hwave_90[t] = np.nanpercentile(Hwave0,90)

# lat_rho = ds['lat_rho'][:]
# lon_rho = ds['lon_rho'][:]
# mask_rho = ds['mask_rho'][:]
ocean_time = ds['ocean_time'][:]

var_list = ['Dwave_avg','Dwave_10','Dwave_25','Dwave_75','Dwave_90',
    'Hwave_avg','Hwave_10','Hwave_25','Hwave_75','Hwave_90',
    'ocean_time']

D = dict()
for var in var_list:
    D[var]=locals()[var]

# outfn = home + 'WQ_data/waves_{:04d}.p'.format(tt)
outfn = home + 'WQ_data/waves_'+suffix+'.p'
pickle.dump(D,open(outfn,'wb'))
