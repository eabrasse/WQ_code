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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

home = '/data0/ebrasseale/'

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

nfiles = len(f_list)

# Time steps are inconsistent across files, so first count 'em up
NT = 0
fn = f_list[0]
ds = nc.Dataset(dir0+fn)

nt,nz,ny,nx = ds['salt'].shape

shorelon = np.zeros((ny))
shorelat = np.zeros((ny))

lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
mask_rho = ds['mask_rho'][:]
x_ind=np.zeros((ny))

for j in range(ny):
    
    # find the edge of the mask
    mask_diff = np.where(np.diff(mask_rho[j,:]))[0]
    
    #if multiple edges, north of TJRE
    if (len(mask_diff)>1)&(lat_rho[j,0]>32.6):
        #look for the edge closest to the previously identified edge
        x_ind[j] = mask_diff[np.argmin(np.abs(x_ind[j-1]-mask_diff))]
        
    #if multiple edges, south of TJRE
    elif (len(mask_diff)>1)&(lat_rho[j,0]<32.6):
        #do outermost edge
        x_ind[j] = mask_diff[0]
        
    elif len(mask_diff)==1:
        x_ind[j] = mask_diff[0]
        
    elif len(mask_diff)==0:
        x_ind[j] = x_ind[j-1]

    shorelon[j] = lon_rho[j,int(x_ind[j])]
    shorelat[j] = lat_rho[j,int(x_ind[j])]
    # x_ind = x_ind0

x_diff = np.diff(shorelon)
#cutoff before jumping at mouth of San Diego bay
cutoff = np.argmax(np.abs(x_diff))
#find where mouth of TJRE is
TJRE_inds = np.where(np.abs(x_diff[:cutoff])>0.0036)[0]
TJ0 = TJRE_inds[0]
TJ1 = TJRE_inds[-1]+1

ds.close()

fig = plt.figure(figsize=(6,8))
ax5 = fig.gca()
ax5.contour(lon_rho,lat_rho,mask_rho,colors='k',levels=[1],linewidths=0.5,alpha=1.0)
ax5.scatter(shorelon[:TJ0],shorelat[:TJ0],color='cornflowerblue')
ax5.scatter(shorelon[TJ1:cutoff],shorelat[TJ1:cutoff],color='cornflowerblue')
yl = ax5.get_ylim()
yav = (yl[0] + yl[1])/2
ax5.set_aspect(1/np.sin(np.pi*yav/180))
ax5.set_ylabel('Latitude')
ax5.set_xlabel('Longitude')

outfn = home + 'WQ_data/shoreline.png'
plt.savefig(outfn)
plt.close('all')
