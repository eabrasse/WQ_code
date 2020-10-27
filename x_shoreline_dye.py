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

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=False
if testing:
    f_list = f_list[:3]

nfiles = len(f_list)

buoylon = -117.169
buoylat = 32.570

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    if NT==0:
        nt,nz,ny,nx = ds['u'].shape
        
        x_list = np.zeros((ny))
        shorelon = np.zeros((ny))
        shorelat = np.zeros((ny))
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        
        refgrid = np.abs(lon_rho-buoylon)+np.abs(lat_rho-buoylat)
        buoy_x = np.where(refgrid==refgrid.min())[1][0]
        buoy_y = np.where(refgrid==refgrid.min())[0][0]

    else:
        nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

dye_01 = np.zeros((NT,ny))
dye_02 = np.zeros((NT,ny))
Dwave = np.zeros((NT))
Hwave = np.zeros((NT))
Lwave = np.zeros((NT))
h = np.zeros((NT))
zeta = np.zeros((NT))
ot = np.zeros((NT))

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # select wave direction and significant wave height
    wetdry_mask_rho = ds['wetdry_mash_rho'][:]
    dye_01_0 = ds['dye_01'][:]
    dye_02_0 = ds['dye_02'][:]
    Dwave0 = ds['Dwave'][:]
    Hwave0 = ds['Hwave'][:]
    Lwave0 = ds['Lwave'][:]
    h0 = ds['h'][:]
    zeta0 = ds['zeta'][:]

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]

    for t in range(nt):
        for j in range(ny):
            x_ind = np.where(wetdry_mask_rho[t,j,:]==0)[0][0]-1
            dye_01[old_nt+t,j] = np.mean(dye_01_0[t,:,j,int(x_ind)],axis=1)
            dye_02[old_nt+t,j] = np.mean(dye_02_0[t,:,j,int(x_ind)],axis=1)
        
    Dwave[old_nt:old_nt+nt] = Dwave0[:,buoy_y,buoy_x]
    Hwave[old_nt:old_nt+nt] = Hwave0[:,buoy_y,buoy_x]
    Lwave[old_nt:old_nt+nt] = Lwave0[:,buoy_y,buoy_x]
    h[old_nt:old_nt+nt] = h0[:,buoy_y,buoy_x]
    zeta[old_nt:old_nt+nt] = zeta0[:,buoy_y,buoy_x]
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1

var_list = ['dye_01','dye_02','Dwave','Hwave','Lwave','ot','lon_rho','lat_rho','mask_rho','h','zeta']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_dye.p'
pickle.dump(D,open(outfn,'wb'))
