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

home = '/dataSIO/ebrasseale/'

dir0 = '/dataSIO/NADB2017/NADB2017_0/Output/'
dir1 = '/dataSIO/NADB2018/'
dir2 = '/dataSIO/NADB2019/'

riv_fn0 = dir0+'river_tracer_4river_NADB2017_0.nc'
riv_fn1 = dir1+'Input/river_tracer_4river_NADB2018.nc'
riv_fn2 = dir2+'Input/river_tracer_4river_NADB2019_0.nc'

rt = np.array([])
Q = np.array([])
for riv_fn in [riv_fn0,riv_fn1,riv_fn2]:
    dsr = nc.Dataset(riv_fn)
    
    rt0 = dsr['river_time'][:]
    rt0 = rt0 * 24*60*60 #river time is in days; match to ocean time in seconds
    rt = np.append(rt,rt0,axis=0)
    
    Q0 = np.sum(dsr['river_transport'][:,:5],axis=1) # indeces 0–4 are TJRE, but they are all zero at the same times
    Q = np.append(Q,Q0,axis=0)
    
    dsr.close()

f_list = []
for my_dir in [dir0,dir1,dir2]:
    f_list0 = os.listdir(my_dir)
    f_list0.sort()
    f_list0 = [my_dir+x for x in f_list0 if x[:9]=='ocean_his']
    f_list.extend(f_list0)

testing=False
if testing:
    f_list = f_list[:3]

nfiles = len(f_list)

buoylon = -117.169
buoylat = 32.570
shorenormal = 265


# Time steps are inconsistent across files, so first count 'em up
NT = {}
NT['total'] = 0
NT['N'] = 0
NT['S'] = 0
NT['NQ'] = 0
NT['SQ'] = 0

for fn in f_list:
    
    ds = nc.Dataset(fn)
    
    if NT['total']==0:
        
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
        nt_NQ = np.sum((Dwave>shorenormal)&(Q[rtt0:rtt1]<-5.0))
        nt_SQ = np.sum((Dwave<shorenormal)&(Q[rtt0:rtt1]<-5.0))

    else:
        nt = ds['ocean_time'].shape[0]
        Dwave0 = ds['Dwave'][:]
        Dwave = Dwave0[:,buoy_y,buoy_x]
        nt_N = np.sum(Dwave>shorenormal)
        nt_S = np.sum(Dwave<shorenormal)
        
        rtt0 = np.argmin(np.abs(rt-ds['ocean_time'][0]))
        rtt1 = np.argmin(np.abs(rt-ds['ocean_time'][-1]))+1
        nt_NQ = np.sum((Dwave>shorenormal)&(Q[rtt0:rtt1]<-5.0))
        nt_SQ = np.sum((Dwave<shorenormal)&(Q[rtt0:rtt1]<-5.0))
        
    NT['total'] += nt
    NT['N'] += nt_N
    NT['S'] += nt_S
    NT['NQ'] += nt_NQ
    NT['SQ'] += nt_SQ
    ds.close()

# surf_dye_01 = {}
# surf_dye_01['N'] = np.zeros((ny,nx))
# surf_dye_01['S'] = np.zeros((ny,nx))
#
# surf_dye_02 = {}
# surf_dye_02['N'] = np.zeros((ny,nx))
# surf_dye_02['S'] = np.zeros((ny,nx))
# surf_dye_02['NQ'] = np.zeros((ny,nx))
# surf_dye_02['SQ'] = np.zeros((ny,nx))
#
# # Now do the extraction and processing
# tt=0
# old_nt = 0
# for fn in f_list:
#     print('file {:d} of {:d}'.format(tt,nfiles))
#     ds = nc.Dataset(dir0+fn)
#
#     # read in dye and wave direction
#
#     dye_01_0 = ds['dye_01'][:]
#     dye_02_0 = ds['dye_02'][:]
#     Dwave0 = ds['Dwave'][:]
#     ot = ds['ocean_time'][:]
#
#     # only consider wave direction at buoy
#     Dwave = Dwave0[:,buoy_y,buoy_x]
#     nt = Dwave.shape[0]
#
#     for t in range(nt):
#
#         rtt = np.argmin(np.abs(rt-ot[t]))
#         # if waves from north...
#         if Dwave[t]>shorenormal:
#             dye_01_N += dye_01_0[t,-1,:,:]/NT_N
#             dye_02_N += dye_02_0[t,-1,:,:]/NT_N
#
#             if Q[rtt]<-5.0:
#                 dye_02_NQ += dye_02_0[t,-1,:,:]/NT_NQ
#
#         # if waves from south...
#         elif Dwave[t]<shorenormal:
#             dye_01_S += dye_01_0[t,-1,:,:]/NT_S
#             dye_02_S += dye_02_0[t,-1,:,:]/NT_S
#
#             if Q[rtt]<-5.0:
#                 dye_02_SQ += dye_02_0[t,-1,:,:]/NT_SQ
#
#
#     ds.close()
#     tt+=1
#
# var_list = ['dye_01_N','dye_01_S','dye_02_N','dye_02_S','dye_02_NQ','dye_02_SQ','NT_N','NT_S','NT_NQ','NT_SQ','NT','lon_rho','lat_rho','mask_rho']
#
# D = dict()
# for var in var_list:
#     D[var]=locals()[var]
#
# outfn = home + 'WQ_data/surface_dye_Q-5.p'
# pickle.dump(D,open(outfn,'wb'))
