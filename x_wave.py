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

home = '/home/ebrasseale/'

dir0 = '/home/x1wu/SDTJRE_EPA/LV4_RUNFILES/Run2017/Couple/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:24]=='ocean_his_LV4_EPA_winter']
nfiles = len(f_list)

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

Dwave_avg = np.zeros((NT))
# Dwave_10 = np.zeros((NT))
# Dwave_25 = np.zeros((NT))
# Dwave_75 = np.zeros((NT))
# Dwave_90 = np.zeros((NT))
Hwave_avg = np.zeros((NT))
# Hwave_10 = np.zeros((NT))
# Hwave_25 = np.zeros((NT))
# Hwave_75 = np.zeros((NT))
# Hwave_90 = np.zeros((NT))
ot = np.zeros((NT))

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # select wave direction and significant wave height

    Dwave = ds['Dwave'][:]
    Hwave = ds['Hwave'][:]
    nt,ny,nx = Dwave.shape

    for t in range(nt):
        # significant wave height is only nonzero along the grid edges
        # find where wave height isn't zero
        y_inds,x_inds = np.where(Hwave[t,:,:]!=0)
        
        Hwave0 = Hwave[t,y_inds,x_inds]
        # since np.nanpercentile doens't handle masked arrays,
        # replace the mask with nan's
        Hwave0 = np.ma.filled(Hwave0,np.nan)
        #repeat for wave direction
        Dwave0 = Dwave[t,y_inds,x_inds]
        Dwave0 = np.ma.filled(Dwave0,np.nan)
        
        # Take average of the remaining indeces
        Dwave_avg[old_nt+t] = np.nanmean(Dwave0)
        # Dwave_10[old_nt+t] = np.nanpercentile(Dwave0,10)
        # Dwave_25[old_nt+t] = np.nanpercentile(Dwave0,25)
        # Dwave_75[old_nt+t] = np.nanpercentile(Dwave0,75)
        # Dwave_90[old_nt+t] = np.nanpercentile(Dwave0,90)
        Hwave_avg[old_nt+t] = np.nanmean(Hwave0)
        # Hwave_10[old_nt+t] = np.nanpercentile(Hwave0,10)
        # Hwave_25[old_nt+t] = np.nanpercentile(Hwave0,25)
        # Hwave_75[old_nt+t] = np.nanpercentile(Hwave0,75)
        # Hwave_90[old_nt+t] = np.nanpercentile(Hwave0,90)
        
    ot[old_nt:old_nt+nt] = ds['ocean_time'][:]
    old_nt += nt
    
    ds.close()
    tt+=1

var_list = ['Dwave_avg','Hwave_avg','ot']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/waves_allyear.p'
pickle.dump(D,open(outfn,'wb'))
