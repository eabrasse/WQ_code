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

buoylon = -117.169
buoylat = 32.570
shorenormal = 259.03

# Time steps are inconsistent across files, so first count 'em up
NT = 0
NT_N = 0
NT_S = 0

for fn in f_list:
    
    ds = nc.Dataset(dir0+fn)
    
    if NT==0:
        nt,nz,ny,nx = ds['salt'].shape
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        
        refgrid = np.abs(lon_rho-buoylon)+np.abs(lat_rho-buoylat)
        buoy_x = np.where(refgrid==refgrid.min())[1][0]
        buoy_y = np.where(refgrid==refgrid.min())[0][0]
        
        Dwave0 = ds['Dwave'][:]
        Dwave = Dwave0[:,buoy_y,buoy_x]
        nt_N = np.sum(Dwave>shorenormal)
        nt_S = np.sum(Dwave<shorenormal)

    else:
        nt = ds['ocean_time'].shape[0]
        Dwave0 = ds['Dwave'][:]
        Dwave = Dwave0[:,buoy_y,buoy_x]
        nt_N = np.sum(Dwave>shorenormal)
        nt_S = np.sum(Dwave<shorenormal)
        
    NT += nt
    NT_N += nt_N
    NT_S += nt_S
    ds.close()

dye_01_N = np.zeros((ny,nx))
dye_01_S = np.zeros((ny,nx))
dye_02_N = np.zeros((ny,nx))
dye_02_S = np.zeros((ny,nx))

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # read in dye and wave direction

    dye_01_0 = ds['dye_01'][:]
    dye_02_0 = ds['dye_02'][:]
    Dwave0 = ds['Dwave'][:]
    
    # only consider wave direction at buoy   
    Dwave = Dwave0[:,buoy_y,buoy_x]
    nt = Dwave.shape[0]
    
    for t in range(nt):
        
        # if waves from north...
        if Dwave[t]>shorenormal:
            dye_01_N += dye_01_0[t,-1,:,:]/NT_N
            dye_02_N += dye_02_0[t,-1,:,:]/NT_N
            
        # if waves from south...
        elif Dwave[t]<shorenormal:
            dye_01_S += dye_01_0[t,-1,:,:]/NT_S
            dye_02_S += dye_02_0[t,-1,:,:]/NT_S 
    
    ds.close()
    tt+=1

var_list = ['dye_01_N','dye_01_S','dye_02_N','dye_02_S','NT_N','NT_S','NT','lon_rho','lat_rho','mask_rho']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/surface_dye.p'
pickle.dump(D,open(outfn,'wb'))
