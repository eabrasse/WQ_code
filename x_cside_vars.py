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
import wqfun

home = '/data0/ebrasseale/'

ind_fn = home+'WQ_data/shore_buoy_inds.p'
Dinds = pickle.load(open(ind_fn,'rb'))
for var in Dinds.keys():
    locals()[var]=Dinds[var]

# year = '2017'
# if year=='2017':
dir2017 = '/data0/NADB2017/NADB2017_0_NEW/'
f_list2017 = os.listdir(dir2017)
f_list2017.sort()
f_list2017 = [x for x in f_list2017 if x[:17]=='ocean_his_NADB_0_']
# elif year=='2018':
dir2018 = '/data0/NADB2018/'
f_list2018 = os.listdir(dir2018)
f_list2018.sort()
f_list2018 = [x for x in f_list2018 if x[:19]=='ocean_his_NADB2018_']

f_list = []
for fn in f_list2017:
    f_list.append(dir2017+fn)
for fn in f_list2018:
    f_list.append(dir2018+fn)


testing=False
if testing:
    f_list = f_list[:1]

ref_depth = 5

nfiles = len(f_list)

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(fn)
    if NT==0:
        # nt,nz,ny,nx = ds['salt'].shape
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        mask_diff = np.zeros((len(jjs)))
        lonshore = np.zeros((len(jjs)))
        latshore = np.zeros((len(jjs)))
        for j in range(len(jjs)):
            mask_diff = np.where(np.diff(mask_rho[jjs[j],:]))[0]
            lonshore[j] = lon_rho[jjs[j],iis[j]]
            latshore[j] = lat_rho[jjs[j],iis[j]]
            
        xshore, yshore = wqfun.ll2xy(lonshore,latshore,lon_rho.min(),lat_rho.min())
        x_rho,y_rho = wqfun.ll2xy(lon_rho,lat_rho,lon_rho.min(),lat_rho.min())
        dxs = np.diff(xshore)
        dys = np.diff(yshore)
        drs = np.sqrt(dxs**2 + dys**2)
        rshore = np.cumsum(drs)
        rshore = np.insert(rshore,0,0,axis=0)


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
u0 = np.zeros((NT,nj))
v0 = np.zeros((NT,nj))
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
    u00 = ds['u'][:]
    v00 = ds['v'][:]
    
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
            u0[old_nt+t,j] = np.nanmean(u00[t,:,jjs[j],int(x_5m_ind):int(x_wd_ind)])
            v0[old_nt+t,j] = np.nanmean(v00[t,:,jjs[j],int(x_5m_ind):int(x_wd_ind)])
        
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


var_list = ['lon_rho','lat_rho','mask_rho','lonshore','latshore',\
'x_rho','y_rho','xshore','yshore','rshore',\
'iis','jjs','iib','jjb','shoreangle',\
'Dwave','Hwave','Lwave','zeta',\
'dye_01','dye_02','u0','v0','ot']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_variables_2017–2018.p'
pickle.dump(D,open(outfn,'wb'))
