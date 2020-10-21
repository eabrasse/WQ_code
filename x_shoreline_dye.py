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

testing=True
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
        
        for j in range(ny):
            x_list[j] = np.where(mask_rho[j,:]==0)[0][0]-2
            shorelon[j] = lon_rho[j,int(x_list[j])]
            shorelat[j] = lat_rho[j,int(x_list[j])]
    else:
        nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

dye_01 = np.zeros((NT,ny))
dye_02 = np.zeros((NT,ny))
Dwave = np.zeros((NT))
Hwave = np.zeros((NT))
ot = np.zeros((NT))

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # select wave direction and significant wave height

    dye_01_0 = ds['dye_01'][:]
    dye_02_0 = ds['dye_02'][:]
    Dwave0 = ds['Dwave'][:]
    Hwave0 = ds['Hwave'][:]

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]

    for j in range(ny):
        dye_01[old_nt:old_nt+nt,j] = np.mean(dye_01_0[:,:,j,int(x_list[j])],axis=1)
        dye_02[old_nt:old_nt+nt,j] = np.mean(dye_02_0[:,:,j,int(x_list[j])],axis=1)
        
    Dwave[old_nt:old_nt+nt] = Dwave0[:,buoy_y,buoy_x]
    Hwave[old_nt:old_nt+nt] = Hwave0[:,buoy_y,buoy_x]
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1

var_list = ['shorelat','shorelon','dye_01','dye_02','Dwave','Hwave','ot','lon_rho','lat_rho','mask_rho']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_dye_-2.p'
pickle.dump(D,open(outfn,'wb'))
