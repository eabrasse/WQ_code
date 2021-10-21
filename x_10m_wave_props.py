#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup

import pickle
import numpy as np
import netCDF4 as nc


ll_fn = '/data0/ebrasseale/WQ_data/shoreline_variables_2017.p'
D=pickle.load(open(ll_fn,'rb'))
var_list = ['x_rho','y_rho','mask_rho','yshore','xshore','rshore','x10m','y10m','h','rshore','lon_rho','lat_rho','iis','iib','jjs','jjb']
for var in var_list:
    locals()[var] = D[var]

fname = '/data0/NADB2017/NADB2017_0_NEW/ocean_his_NADB_0_new_00001.nc'
ds = nc.Dataset(fname)
nt = ds['ocean_time'][:].shape
Hwave = ds['Hwave'][:]
Dwave = ds['Dwave'][:]
zeta = ds['zeta']

r_refs = [9000,12500,17000]
D = {}
for ref in r_refs:
    D[ref] = {}
    j = np.argmin(np.abs(rshore-ref))
    js = [jjb[j],jjs[j]]
    jj0 = int(np.min(js))
    jj1 = int(np.max(js))+1
    jvec = range(jj0,jj1)
    ni = iis[j]-iib[j]
    
    D[ref]['rdist'] = np.zeros((ni))
    D[ref]['Hwave_vec'] = np.zeros((nt,ni))
    D[ref]['Dwave_vec'] = np.zeros((nt,ni))
    D[ref]['zeta_vec'] = np.zeros((nt,ni))
    #look at each x index between the buoy and the shore
    for ii in range(ni):
        xvec = x_rho[jvec,iib[j]+ii]
        yvec = y_rho[jvec,iib[j]+ii]
        xdiff = xvec-xshore[j]
        ydiff = yvec-yshore[j]
        
        local_angle = np.arctan2(xdiff,ydiff)*180/np.pi+90
        angle_diff = np.abs(local_angle-shoreangle[j])
        jjind = np.argmin(angle_diff)
        
        # hvec[ii] = -h[jvec[jjind],iib[j]+ii]
        D[ref]['rdist'][ii] = np.sqrt((x_rho[jvec[jjind],iib[j]+ii]-xshore[j])**2+(y_rho[jvec[jjind],iib[j]+ii]-yshore[j])**2)

        D[ref]['Hwave_vec'][:,ii] = Hwave[:,jvec[jjind],iib[j]+ii]
        D[ref]['Dwave_vec'][:,ii] = Dwave[:,jvec[jjind],iib[j]+ii]
        D[ref]['zeta_vec'][:,ii] = zeta[:,jvec[jjind],iib[j]+ii]

outfn = '/data0/ebrasseale/WQ_data/10m_wave_props.p'
pickle.dump(D,outfn)

ds.close()