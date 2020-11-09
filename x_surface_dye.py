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

riv_fn = dir0+'river_tracer_4river_NADB2017_0.nc'
dsr = nc.Dataset(riv_fn)
rt = dsr['river_time'][:]
rt = rt * 24*60*60 #river time is in days; match to ocean time in seconds
Q = dsr['river_transport'][:,0] # indeces 0–4 are TJRE, but they are all zero at the same times

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
NT_NQ = 0
NT_SQ = 0

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
        
        rtt0 = np.argmin(np.abs(rt-ds['ocean_time'][0]))
        rtt1 = np.argmin(np.abs(rt-ds['ocean_time'][-1]))+1
        nt_NQ = np.sum((Dwave>shorenormal)&(Q[rtt0:rtt1]!=0))
        nt_SQ = np.sum((Dwave<shorenormal)&(Q[rtt0:rtt1]!=0))

    else:
        nt = ds['ocean_time'].shape[0]
        Dwave0 = ds['Dwave'][:]
        Dwave = Dwave0[:,buoy_y,buoy_x]
        nt_N = np.sum(Dwave>shorenormal)
        nt_S = np.sum(Dwave<shorenormal)
        
        rtt0 = np.argmin(np.abs(rt-ds['ocean_time'][0]))
        rtt1 = np.argmin(np.abs(rt-ds['ocean_time'][-1]))+1
        nt_NQ = np.sum((Dwave>shorenormal)&(Q[rtt0:rtt1]!=0))
        nt_SQ = np.sum((Dwave<shorenormal)&(Q[rtt0:rtt1]!=0))
        
    NT += nt
    NT_N += nt_N
    NT_S += nt_S
    NT_NQ += nt_NQ
    NT_SQ += nt_SQ
    ds.close()

dye_01_N = np.zeros((ny,nx))
dye_01_S = np.zeros((ny,nx))
dye_02_N = np.zeros((ny,nx))
dye_02_S = np.zeros((ny,nx))
dye_02_NQ = np.zeros((ny,nx))
dye_02_SQ = np.zeros((ny,nx))

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
    ot = ds['ocean_time'][:]
    
    # only consider wave direction at buoy   
    Dwave = Dwave0[:,buoy_y,buoy_x]
    nt = Dwave.shape[0]
    
    for t in range(nt):

        rtt = np.argmin(np.abs(rt-ot[t]))
        # if waves from north...
        if Dwave[t]>shorenormal:
            dye_01_N += dye_01_0[t,-1,:,:]/NT_N
            dye_02_N += dye_02_0[t,-1,:,:]/NT_N
            
            if Q[rtt]!=0:
                dye_02_NQ += dye_02_0[t,-1,:,:]/NT_NQ
            
        # if waves from south...
        elif Dwave[t]<shorenormal:
            dye_01_S += dye_01_0[t,-1,:,:]/NT_S
            dye_02_S += dye_02_0[t,-1,:,:]/NT_S      
                  
            if Q[rtt]!=0:
                dye_02_SQ += dye_02_0[t,-1,:,:]/NT_SQ

    
    ds.close()
    tt+=1

var_list = ['dye_01_N','dye_01_S','dye_02_N','dye_02_S','dye_02_NQ','dye_02_SQ','NT_N','NT_S','NT_NQ','NT_SQ','NT','lon_rho','lat_rho','mask_rho']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/surface_dye_Q.p'
pickle.dump(D,open(outfn,'wb'))
