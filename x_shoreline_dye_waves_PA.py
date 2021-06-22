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

ind_fn = home+'WQ_data/shore_buoy_inds_PA.p'
Dinds = pickle.load(open(ind_fn,'rb'))
for var in Dinds.keys():
    locals()[var]=Dinds[var]

dir0 = '/data0/NADB2017/NADB2017_0_NEW/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=False
if testing:
    f_list = f_list[:1]

ref_depth = 5

nfiles = len(f_list)

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(dir0+fn)
    if NT==0:
        # nt,nz,ny,nx = ds['salt'].shape
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        mask_diff = np.zeros((len(jjs)))
        for j in range(len(jjs)):
            mask_diff = np.where(np.diff(mask_rho[jjs[j],:]))[0]


    nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()


nj = len(jjs)

dye_01 = np.zeros((NT,nj))
dye_02 = np.zeros((NT,nj))
Dwave = np.zeros((NT,nj))
Hwave = np.zeros((NT,nj))
Lwave = np.zeros((NT,nj))
zeta = np.zeros((NT,nj))
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
    
    H = h0+zeta0

    ocean_time = ds['ocean_time'][:]
    nt = ocean_time.shape[0]


    for j in range(nj):
        #loop through time steps, 
        # because extraction indexes depend on time-varying wetdry mask
        for t in range(nt):
            # find the edge of the mask
            wd_mask_diff = np.where(np.diff(wetdry_mask_rho[t,jjs[j],:]))[0]
            #find where depth crosses from deeper than ref_depth to shallower
            depth_diff = np.where(np.diff(np.sign(H[t,jjs[j],:]-ref_depth)))[0]
    
            #if multiple edges, north of TJRE
            if (len(mask_diff)>1)&(lat_rho[jjs[j],0]>32.6):
                #look for the edge closest to the previously identified edge
                x_wd_ind = wd_mask_diff[np.argmin(np.abs(x_wd_ind-wd_mask_diff))]
                x_5m_ind = depth_diff[np.argmin(np.abs(x_5m_ind-depth_diff))]

            #if multiple edges, south of TJRE
            elif (len(mask_diff)>1)&(lat_rho[jjs[j],0]<32.6):
                #do outermost edge
                x_wd_ind = wd_mask_diff[0]
                x_5m_ind = depth_diff[0]

            elif len(mask_diff)==1:
                x_wd_ind = wd_mask_diff[0]
                x_5m_ind = depth_diff[0]

            #go offshore of the wet/dry mask by a tad
            # x_wd_ind = x_wd_ind - 2
            dye_01[old_nt+t,j] = np.nanmean(dye_01_0[t,:,jjs[j],int(x_5m_ind):int(x_wd_ind)])
            dye_02[old_nt+t,j] = np.nanmean(dye_02_0[t,:,jjs[j],int(x_5m_ind):int(x_wd_ind)])
        
        #the wave extractions are from the same points in space no matter what the time
        #so they can exist outside of the loop
        Dwave[old_nt:old_nt+nt,j] = Dwave0[:,jjb[j],iib[j]]
        Hwave[old_nt:old_nt+nt,j] = Hwave0[:,jjb[j],iib[j]]
        Lwave[old_nt:old_nt+nt,j] = Lwave0[:,jjb[j],iib[j]]
        zeta[old_nt:old_nt+nt,j] = zeta0[:,jjb[j],iib[j]]
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1


var_list = ['dye_01','dye_02','ot','lon_rho','lat_rho','mask_rho','iis','jjs','iib','jjb','Dwave','Hwave','Lwave','zeta','shoreangle']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_dye_waves_05m_PA.p'
pickle.dump(D,open(outfn,'wb'))
