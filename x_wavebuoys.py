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

testing=False
if testing:
    f_list = f_list[:3]

nfiles = len(f_list)

# buoyname = ['TJRE','PB']
# buoylon = [-117.169,-117.15]
# buoylat = [32.570,32.446]
href = 10

# nbuoys = len(buoylon)
nbins=10
nbuoys=nbins-1
buoy_x = np.zeros((nbuoys))
buoy_y = np.zeros((nbuoys))
buoylon = np.zeros((nbuoys))
buoylat = np.zeros((nbuoys))
h = np.zeros((nbuoys))

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    nt = ds['ocean_time'].shape[0]
    if NT==0:
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        
        ny = lat_rho.shape[0]
        ss = int(np.floor(ny/nbins))
        ss2 = int(np.floor(ss/2))
        
        for b in range(nbuoys):
            buoy_y[b] = b*ss+ss2
            
            hvec = h0[b*ss+ss2,:]
            buoy_x[b] = np.argmin(np.abs(hvec-href))
            j = int(buoy_y[b])
            i = int(buoy_x[b])
            buoylon[b]=lon_rho[j,i]
            buoylat[b]=lat_rho[j,i]
            h[b] = h0[j,i]
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
    zeta0 = ds['zeta'][:]

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

outfn = home + 'WQ_data/wavebuoys_{:d}buoys.p'.format(nbuoys)
pickle.dump(D,open(outfn,'wb'))
