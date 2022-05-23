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
import time

home = '/dataSIO/ebrasseale/'

dir0 = '/dataSIO/NADB2017/NADB2017_0_NEW/'
dir1 = '/dataSIO/NADB2018/'
dir2 = '/dataSIO/NADB2019/'

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
nt_guess = 240*nfiles
ot = np.zeros((nt_guess))

# build variable list to filter
ds = nc.Dataset(f_list[0])
var_2gf_list = [] #variables to godin filter
var_ngf_list = [] #variables to not godin filter
for var in ds.variables.keys():
    var0 = ds[var]
    try:
        if 'ocean_time' in var0.coordinates:
            var_2gf_list.append(var) # if time-varying, add name to "to godin filter" list
        else:
            var_ngf_list.append(var) # if not time-varying, add name to "not godin filter" list
    except:
        var_ngf_list.append(var) #variable doesn't have 'coordinates'. save it as "not godin filter"

#
# tic = time.perf_counter()
# print('Counting up time steps and reading in surface dye...')
# count =0
# for fn in f_list:
#
#     printout = f'(working on file {count} of {nfiles})'
#     print(printout)
#
#     ds = nc.Dataset(fn)
#
#     if count==0:
#
#         nt,nz,ny,nx = ds['salt'].shape
#
#         lon_rho = ds['lon_rho'][:]
#         lat_rho = ds['lat_rho'][:]
#         mask_rho = ds['mask_rho'][:]
#
#         refgrid = np.abs(lon_rho-buoylon)+np.abs(lat_rho-buoylat)
#         buoy_x = np.where(refgrid==refgrid.min())[1][0]
#         buoy_y = np.where(refgrid==refgrid.min())[0][0]
#
#         surf_dye_01 = {}
#         surf_dye_02 = {}
#         for cond in cond_list:
#             surf_dye_01[cond] = np.zeros((ny,nx))
#             surf_dye_02[cond] = np.zeros((ny,nx))
#
#     else:
#         nt = ds['ocean_time'].shape[0]
#
#     # count time steps for each condition
#     Dwave0 = ds['Dwave'][:]
#     Dwave = Dwave0[:,buoy_y,buoy_x]
#     nt_N = np.sum(Dwave>shorenormal)
#     nt_S = np.sum(Dwave<shorenormal)
#
#     rtt0 = np.argmin(np.abs(rt-ds['ocean_time'][0]))
#     rtt1 = np.argmin(np.abs(rt-ds['ocean_time'][-1]))+1
#     nt_NQ = np.sum((Dwave>shorenormal)&(Q[rtt0:rtt1]<-5.0))
#     nt_SQ = np.sum((Dwave<shorenormal)&(Q[rtt0:rtt1]<-5.0))
#
#     ot0 = ds['ocean_time'][:]
#     ot[NT['total']:NT['total']+nt] = ot0
#     NT['total'] += nt
#     NT['N'] += nt_N
#     NT['S'] += nt_S
#     NT['NQ'] += nt_NQ
#     NT['SQ'] += nt_SQ
#
#     # add up sea surface
#     dye_01_0 = ds['dye_01'][:]
#     dye_02_0 = ds['dye_02'][:]
#
#     for t in range(nt):
#
#         rtt = np.argmin(np.abs(rt-ot0[t]))
#         # if waves from north...
#         if Dwave[t]>shorenormal:
#             surf_dye_01['N'] += dye_01_0[t,-1,:,:]
#             surf_dye_02['N'] += dye_02_0[t,-1,:,:]
#
#             if Q[rtt]<-5.0:
#                 surf_dye_02['NQ'] += dye_02_0[t,-1,:,:]
#
#         # if waves from south...
#         elif Dwave[t]<shorenormal:
#             surf_dye_01['S'] += dye_01_0[t,-1,:,:]
#             surf_dye_02['S'] += dye_02_0[t,-1,:,:]
#
#             if Q[rtt]<-5.0:
#                 surf_dye_02['SQ'] += dye_02_0[t,-1,:,:]
#
#
#     ds.close()
#     count+=1
#
# print('Counting time steps and reading surface dye complete!')
# toc = time.perf_counter()
# count_time = f"Counting time steps took {toc-tic:0.4f} seconds"
# print(count_time)
#
# tic = time.perf_counter()
# print('Dividing cumulative surface dye by time steps to get mean...')
# # divide by time steps to get mean
# for cond in cond_list:
#     surf_dye_01[cond] = surf_dye_01[cond]/NT[cond]
#     surf_dye_02[cond] = surf_dye_02[cond]/NT[cond]
# toc = time.perf_counter()
# count_time = f"Dividing dye took {toc-tic:0.4f} seconds"
# print(count_time)
#
# ot = ot[ot>0]
# var_list = ['surf_dye_01','surf_dye_02','NT','Dwave','Q','ot','rt','lon_rho','lat_rho','mask_rho']
#
# D = dict()
# for var in var_list:
#     D[var]=locals()[var]
#
# outfn = home + 'WQ_data/surface_dye_2017-2019_Q-5.p'
# print(f'Saving to {outfn}')
# pickle.dump(D,open(outfn,'wb'))
