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
        nt,nz,ny,nx = ds['salt'].shape
        
        shorelon = np.zeros((ny))
        shorelat = np.zeros((ny))
        x_ind=np.zeros((ny))
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        
        refgrid = np.abs(lon_rho-buoylon)+np.abs(lat_rho-buoylat)
        buoy_x = np.where(refgrid==refgrid.min())[1][0]
        buoy_y = np.where(refgrid==refgrid.min())[0][0]
        h = h0[buoy_y,buoy_x]
        
        # for j in range(ny):
        #     x_ind = np.where(mask_rho[j,:]==0)[0][0]-1
        #     if lon_rho[j,int(x_ind)]<-117.2:
        #         x_ind = np.where(mask_rho[j,:])[0][2]-1
        #     shorelon[j] = lon_rho[j,int(x_ind)]
        #     shorelat[j] = lat_rho[j,int(x_ind)]
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
    wetdry_mask_rho = ds['wetdry_mask_rho'][:]
    dye_01_0 = ds['dye_01'][:]
    dye_02_0 = ds['dye_02'][:]
    Dwave0 = ds['Dwave'][:]
    Hwave0 = ds['Hwave'][:]
    Lwave0 = ds['Lwave'][:]
    zeta0 = ds['zeta'][:]

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]

    for t in range(nt):
        for j in range(ny):
            # find the edge of the mask
            mask_diff = np.where(np.diff(wetdry_mask_rho[t,j,:]))[0]
    
            #if multiple edges, north of TJRE
            if (len(mask_diff)>1)&(lat_rho[j,0]>32.6):
                #look for the edge closest to the previously identified edge
                x_wd_ind = mask_diff[np.argmin(np.abs(x_wd_ind-mask_diff))]
        
            #if multiple edges, south of TJRE
            elif (len(mask_diff)>1)&(lat_rho[j,0]<32.6):
                #do outermost edge
                x_wd_ind = mask_diff[0]
        
            elif len(mask_diff)==1:
                x_wd_ind = mask_diff[0]

            dye_01[old_nt+t,j] = np.mean(dye_01_0[t,:,j,int(x_wd_ind)])
            dye_02[old_nt+t,j] = np.mean(dye_02_0[t,:,j,int(x_wd_ind)])
        
    Dwave[old_nt:old_nt+nt] = Dwave0[:,buoy_y,buoy_x]
    Hwave[old_nt:old_nt+nt] = Hwave0[:,buoy_y,buoy_x]
    Lwave[old_nt:old_nt+nt] = Lwave0[:,buoy_y,buoy_x]
    zeta[old_nt:old_nt+nt] = zeta0[:,buoy_y,buoy_x]
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1

#trim the top of the data, since it doesn't follow the shoreline
x_diff = np.diff(shorelon)
cutoff = np.argmax(np.abs(x_diff))
shorelon = shorelon[:cutoff]
shorelat = shorelat[:cutoff]
dye_01 = dye_01[:,:cutoff]
dye_02 = dye_02[:,:cutoff]

var_list = ['dye_01','dye_02','Dwave','Hwave','Lwave','ot','lon_rho','lat_rho','mask_rho','h','zeta','shorelon','shorelat']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_dye.p'
pickle.dump(D,open(outfn,'wb'))
