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

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=True
if testing:
    f_list = f_list[:3]

nfiles = len(f_list)

buoylon = [-117.169,-117.15]
buoylat = [32.570,32.446]

nbuoys = len(buoylon)
buoy_x = np.zeros((nbuoys))
buoy_y = np.zeros((nbuoys))
h = np.zeros((nbuoys))

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    if NT==0:
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        
        for b in range(nbuoys):
            refgrid = np.abs(lon_rho-buoylon[b])+np.abs(lat_rho-buoylat[b])
            buoy_x[b] = np.where(refgrid==refgrid.min())[1][0]
            buoy_y[b] = np.where(refgrid==refgrid.min())[0][0]
            h[b] = h0[int(buoy_y[b]),int(buoy_x[b])]

    else:
        nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

Dwave = np.zeros((nbuoys,NT))
Hwave = np.zeros((nbuoys,NT))
Lwave = np.zeros((nbuoys,NT))
zeta = np.zeros((nbuoys,NT))
ot = np.zeros((NT))

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # select wave direction and significant wave height
    Dwave0 = ds['Dwave'][:]
    Hwave0 = ds['Hwave'][:]
    Lwave0 = ds['Lwave'][:]
    

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]

    for b in range(nbuoys):
        j = int(buoy_y[b])
        i = int(buoy_x[b])
        Dwave[b,old_nt:old_nt+nt] = Dwave0[:,j,i]
        Hwave[b,old_nt:old_nt+nt] = Hwave0[:,j,i]
        Lwave[b,old_nt:old_nt+nt] = Lwave0[:,j,i]
        zeta[b,old_nt:old_nt+nt] = zeta0[:,j,i]
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1

var_list = ['buoylon','buoylat','Dwave','Hwave','Lwave','ot','lon_rho','lat_rho','mask_rho','h','zeta']
#adding a note here to commit
D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/wavebuoys_TJRE_PB.p'
pickle.dump(D,open(outfn,'wb'))
