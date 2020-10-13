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
f_list = [x for x in f_list if x[:5]=='ocean']
nfiles = len(f_list)

tt=0
for fn in f_list:
    print('day {} of {}'.format(tt,nfiles))
    ds = nc.Dataset(fn)
    
    if tt==0:
        #initialize time series
        nt = ds['ocean_time'].shape()[0]
        
        NT = nt*nfiles

        Dwave_avg = np.zeros((NT))
        Dwave_10 = np.zeros((NT))
        Dwave_25 = np.zeros((NT))
        Dwave_75 = np.zeros((NT))
        Dwave_90 = np.zeros((NT))
        Hwave_avg = np.zeros((NT))
        Hwave_10 = np.zeros((NT))
        Hwave_25 = np.zeros((NT))
        Hwave_75 = np.zeros((NT))
        Hwave_90 = np.zeros((NT))
        ot = np.zeros((NT))


    # select wave direction and significant wave height

    Dwave = ds['Dwave'][:]
    Hwave = ds['Hwave'][:]

    for t in range(nt):
        y_inds,x_inds = np.where(Hwave[t,:,:]!=0)
        Hwave0 = Hwave[t,y_inds,x_inds]
        Hwave0 = np.ma.filled(Hwave0,np.nan)
        Dwave0 = Dwave[t,y_inds,x_inds]
        Dwave0 = np.ma.filled(Dwave0,np.nan)
    
        Dwave_avg[tt*nt+t] = np.nanmean(Dwave0)
        Dwave_10[tt*nt+t] = np.nanpercentile(Dwave0,10)
        Dwave_25[tt*nt+t] = np.nanpercentile(Dwave0,25)
        Dwave_75[tt*nt+t] = np.nanpercentile(Dwave0,75)
        Dwave_90[tt*nt+t] = np.nanpercentile(Dwave0,90)
        Hwave_avg[tt*nt+t] = np.nanmean(Hwave0)
        Hwave_10[tt*nt+t] = np.nanpercentile(Hwave0,10)
        Hwave_25[tt*nt+t] = np.nanpercentile(Hwave0,25)
        Hwave_75[tt*nt+t] = np.nanpercentile(Hwave0,75)
        Hwave_90[tt*nt+t] = np.nanpercentile(Hwave0,90)

    ot[tt*nt:(tt+1)*nt] = ds['ocean_time'][:]
    
    ds.close()
    tt+=1

var_list = ['Dwave_avg','Dwave_10','Dwave_25','Dwave_75','Dwave_90',
    'Hwave_avg','Hwave_10','Hwave_25','Hwave_75','Hwave_90',
    'ocean_time']

D = dict()
for var in var_list:
    D[var]=locals()[var]

# outfn = home + 'WQ_data/waves_{:04d}.p'.format(tt)
outfn = home + 'WQ_data/waves_allyear.p'
pickle.dump(D,open(outfn,'wb'))
