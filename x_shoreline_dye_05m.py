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

dir0 = '/data0/NADB2017/NADB2017_0_NEW/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=False
if testing:
    f_list = f_list[:3]

ref_depth = 5

nfiles = len(f_list)
#
# buoylon = -117.169
# buoylat = 32.570

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
        
        # refgrid = np.abs(lon_rho-buoylon)+np.abs(lat_rho-buoylat)
        # buoy_x = np.where(refgrid==refgrid.min())[1][0]
        # buoy_y = np.where(refgrid==refgrid.min())[0][0]
        # h = h0[buoy_y,buoy_x]
        
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

    else:
        nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()

#trim the top of the data, since it doesn't follow the shoreline
#find where shoreline jumps
x_diff = np.diff(shorelon)
#identify largest jump
cutoff = np.argmax(np.abs(x_diff))

#and remove TJRE mouth
TJRE_inds = np.where(np.abs(x_diff[:cutoff])>0.0036)[0]
TJ0 = TJRE_inds[0]
TJ1 = TJRE_inds[-1]+1
j_inds = list(range(TJ0))+list(range(TJ1,cutoff))

#trim shoreline data
shorelon = shorelon[j_inds]
shorelat = shorelat[j_inds]

nj = len(j_inds)

dye_01 = np.zeros((NT,nj))
dye_02 = np.zeros((NT,nj))
# Dwave = np.zeros((NT))
# Hwave = np.zeros((NT))
# Lwave = np.zeros((NT))
# h = np.zeros((NT))
# zeta = np.zeros((NT))
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
    # Dwave0 = ds['Dwave'][:]
    # Hwave0 = ds['Hwave'][:]
    # Lwave0 = ds['Lwave'][:]
    zeta0 = ds['zeta'][:]
    
    H = h0+zeta0

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]

    for t in range(nt):
        j_count=0
        for j in j_inds:
            # find the edge of the mask
            wd_mask_diff = np.where(np.diff(wetdry_mask_rho[t,j,:]))[0]
            #find where depth crosses from deeper than ref_depth to shallower
            depth_diff = np.where(np.diff(np.sign(H[t,j,:]-ref_depth)))[0]
            # x_1m_ind = np.argmin(np.abs(H[t,j,:]-ref_depth))
    
            #if multiple edges, north of TJRE
            if (len(mask_diff)>1)&(lat_rho[j,0]>32.6):
                #look for the edge closest to the previously identified edge
                x_wd_ind = wd_mask_diff[np.argmin(np.abs(x_wd_ind-wd_mask_diff))]
                x_1m_ind = depth_diff[np.argmin(np.abs(x_1m_ind-depth_diff))]

            #if multiple edges, south of TJRE
            elif (len(mask_diff)>1)&(lat_rho[j,0]<32.6):
                #do outermost edge
                x_wd_ind = wd_mask_diff[0]
                x_1m_ind = depth_diff[0]

            elif len(mask_diff)==1:
                x_wd_ind = wd_mask_diff[0]
                x_1m_ind = depth_diff[0]

            #go offshore of the wet/dry mask by a tad
            # x_wd_ind = x_wd_ind - 2
            dye_01[old_nt+t,j_count] = np.nanmean(dye_01_0[t,:,j,int(x_1m_ind):int(x_wd_ind)])
            dye_02[old_nt+t,j_count] = np.nanmean(dye_02_0[t,:,j,int(x_1m_ind):int(x_wd_ind)])
            
            j_count+=1
        
    # Dwave[old_nt:old_nt+nt] = Dwave0[:,buoy_y,buoy_x]
    # Hwave[old_nt:old_nt+nt] = Hwave0[:,buoy_y,buoy_x]
    # Lwave[old_nt:old_nt+nt] = Lwave0[:,buoy_y,buoy_x]
    # zeta[old_nt:old_nt+nt] = zeta0[:,buoy_y,buoy_x]
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1


var_list = ['dye_01','dye_02','ot','lon_rho','lat_rho','mask_rho','shorelon','shorelat']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_dye_depth_05m.p'
pickle.dump(D,open(outfn,'wb'))
