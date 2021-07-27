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

# dir0 = '/data0/NADB2017/NADB2017_0_NEW/'
dir0 = '/home/x1wu/SDTJRE_EPA/mfiles/Run2016_2017/NAM_forcing/LV4/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:13]=='roms_nam_LV4_']
nfiles = len(f_list)

# buoyname = ['TJRE','PB']
# buoylon = [-117.169,-117.15]
# buoylat = [32.570,32.446]
buoylon = -117.169
buoylat = 32.570

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    if NT==0: #first index
        lonr = ds['lon'][:]
        latr = ds['lat'][:]
        latlondiff = np.sqrt((lonr-buoylon)**2 + (latr-buoylat)**2)
        iind = np.argwhere(latlondiff==latlondiff.min())[0][1]
        jind = np.argwhere(latlondiff==latlondiff.min())[0][0]
    nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

Uwind = np.zeros((NT))
Vwind = np.zeros((NT))
ot = np.zeros((NT))



# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(dir0+fn)

    # select wave direction and significant wave height

    Uwind0 = ds['Uwind'][:]
    # Vwind0 = ds['vwind'][:]
    nt,ny,nx = Uwind.shape


    Uwind[old_nt:old_nt+nt] = Uwind0
    Vwind[old_nt:old_nt+nt] = ds['Vwind'][:]
    ot[old_nt:old_nt+nt] = ds['ocean_time'][:]
    old_nt += nt
    
    ds.close()
    tt+=1

var_list = ['Uwind','Vwind','ot']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/wind_near_TJRE_2017.p'
pickle.dump(D,open(outfn,'wb'))
