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

# code to compare Punta Bandera river input file with model behavior
# Qin should = Qout, dye_in should = dye_out, etc.

home = '/data0/ebrasseale/'

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
riv_fn = [x for x in f_list if x[0]=='r'][0]
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=True
if testing:
    f_list = f_list[0]

# nfiles = len(f_list)
dsr = nc.Dataset(dir0+riv_fn)
i0 = dsr['river_Xposition'][5]
j0 = dsr['river_Eposition'][5]

i0 = int(i0)
j0 = int(j0)

Q_in = dsr['river_transport'][:,5] # dimensions: river_time, river
dye_in = dsr['river_dye_01'][:,:,5] # dimensions of: river_time, s_rho, river
# I checked beforehand and saw that river_dye_01 is vertically homogeneous
Qdye_in = Q_in * np.mean(dye_in,axis=1)

# count=0
# for fn in f_list:
fn = f_list[0]
ds = nc.Dataset(dir0+fn)
# if count == 0:
lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
mask_rho = ds['mask_rho'][:]

# derive river transport at point
h = ds['h'][:]
ubar = ds['ubar'][:]
pn = ds['pn'][:]
dy = 1/pn[j0,i0]

ubar_PB = ubar[:,j0,i0] #should have dimensions of time
Q_out = ubar_PB * h[j0,i0] * dy

dye_01 = ds['dye_01'][:]
dye_out = dye_01[:,:,j0,i0] # dimensions of time, z
Qdye_out = Q_out * np.mean(dye_out)
    
# count+=1
river_time = dsr['river_time'][:]
ocean_time = ds['ocean_time'][:]

var_list = ['lon_rho','lat_rho','mask_rho','h','dy','i0','j0',
    'Q_in','dye_in','Qdye_in',
    'Q_out','dye_out','Qdye_out',
    'river_time','ocean_time']

D = dict()
for var in var_list:
    D[var]=locals()[var]
    
ds.close()
dsr.close()

outfn = home + 'WQ_data/x_rivers_in_vs_out.p'
pickle.dump(D,open(outfn,'wb'))
