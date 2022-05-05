#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 Code to test extraction of variables from the COAWST model
 using a variable surf zone width
 it only extracts waves at one location (in line with NDBC buoy Station 46235 - Imperial Beach Nearshore, CA (155))
 to simulate how I expect the model to be adapted for other users

 Note: for manuscript, I ended up using the 5-m isobath and not the surf zone.
 However, this code extracted wave properties at the offshore wave buoy which I ended up using.
 I also only used data through 2017 for analysis and this extracts 2017 and 2018.
"""

# setup
import os
import sys
import pickle
import numpy as np
import netCDF4 as nc
import wqfun



home = '/dataSIO/ebrasseale/'


dir2017 = '/dataSIO/NADB2017/NADB2017_0_NEW/'
dir2018 = '/dataSIO/NADB2018/'
dir2019 = '/dataSIO/NADB2019/'
f_list = []
for my_dir in [dir2017,dir2018,dir2019]:
    f_list0 = os.listdir(my_dir)
    f_list0.sort()
    f_list0 = [my_dir+x for x in f_list0 if x[:9]=='ocean_his']
    f_list.extend(f_list0)



iis,jjs = wqfun.get_shore_inds(f_list[0])

testing=False
if testing:
    f_list = f_list[:1]

ref_depth = 5
buoylat = 32.570
buoylon = -117.169

nfiles = len(f_list)

# Time steps are inconsistent across files, so first count 'em up
NT = 0
for fn in f_list:
    ds = nc.Dataset(fn)
    if NT==0:
        # on the first file, load in some grid data that
        # won't change in the extraction
        
        lon_rho = ds['lon_rho'][:]
        lat_rho = ds['lat_rho'][:]
        mask_rho = ds['mask_rho'][:]
        h0 = ds['h'][:]
        lonshore = np.zeros((len(jjs)))
        latshore = np.zeros((len(jjs)))
        for j in range(len(jjs)):
            lonshore[j] = lon_rho[jjs[j],iis[j]]
            latshore[j] = lat_rho[jjs[j],iis[j]]
            
        xshore, yshore = wqfun.ll2xy(lonshore,latshore,lon_rho.min(),lat_rho.min())
        x_rho,y_rho = wqfun.ll2xy(lon_rho,lat_rho,lon_rho.min(),lat_rho.min())
        dxs = np.diff(xshore)
        dys = np.diff(yshore)
        drs = np.sqrt(dxs**2 + dys**2)
        rshore = np.cumsum(drs)
        rshore = np.insert(rshore,0,0,axis=0)
        
        #calculate buoy index
        bi, bj = wqfun.find_nearest_ind_2D(lon_rho,lat_rho,buoylon,buoylat)


    nt = ds['ocean_time'].shape[0]
    NT += nt
    ds.close()


nj = len(jjs)

# shoreline extracted values
dye_01 = np.zeros((NT,nj)) # PB dye
dye_02 = np.zeros((NT,nj)) # TJRE dye
u0 = np.zeros((NT,nj)) # velocity in one direction (eight degrees from E-W)
v0 = np.zeros((NT,nj)) # velocity in an orthogonal direction (eight degreees from N-S)

# buoy wave values
Dwave = np.zeros((NT)) # mean wave direction
Hwave = np.zeros((NT)) # significant wave height
Lwave = np.zeros((NT)) # mean wave length
zeta = np.zeros((NT)) # SSH

# time
ot = np.zeros((NT)) # time

#prepare index finding dictionaries
WD_rho = {}
WD_rho['ind'] = 0

WD_u = {}
WD_u['ind'] = 0

WD_v = {}
WD_v['ind'] = 0

X5m = {}
X5m['ind'] = 0

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print(f'file {(tt+1):d} of {nfiles:d}')
    ds = nc.Dataset(fn)

    # select wave direction and significant wave height
    WD_rho['data'] = ds['wetdry_mask_rho'][:]
    WD_u['data'] = ds['wetdry_mask_u'][:]
    WD_v['data'] = ds['wetdry_mask_v'][:]
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

    #the wave extractions are from the same points in space no matter what the time
    #so they can exist outside of the shoreline indexing loop
    Dwave[old_nt:old_nt+nt] = Dwave0[:,bj,bi]
    Hwave[old_nt:old_nt+nt] = Hwave0[:,bj,bi]
    Lwave[old_nt:old_nt+nt] = Lwave0[:,bj,bi]
    zeta[old_nt:old_nt+nt] = zeta0[:,bj,bi]

    for j in range(nj):
        print(f'j = {(j+1)} of {nj}')
        #loop through time steps, 
        # because extraction indexes depend on time-varying wetdry mask
        north_of_TJRE = lat_rho[jjs[j],0]>32.6
        for t in range(nt):
            
            for WD in WD_rho,WD_u,WD_v:
                WD['diff'] = np.where(np.diff(WD['data'][t,jjs[j],:(iis[j]+2)])<0)[0]
            
            X5m['diff'] = np.where(np.diff(np.sign(H[t,jjs[j],:(iis[j]+2)]-ref_depth)))[0]
        
        
            for DD in WD_rho,WD_u,WD_v,X5m:
                if len(DD['diff'])==1:
                    DD['ind'] = DD['diff'][0]
                elif len(DD['diff'])>1:
                    # south of TJRE, just choose the most western point
                    # north of TJRE, San Diego Bay and Pt Loma get tricky
                    if north_of_TJRE:
                        # to avoid a weird problem where sometimes the wetdry mask leaves an isolated wet
                        # spot in an otherwise dry tidal zone, make sure we're only considering wet areas with 2 wet indexes west
                        # of them
                        if DD in [WD_rho,WD_u,WD_v]:
                            new_diff_list = []
                            for ind_diff in DD['diff']:
                                if ~np.any(DD['data'][t,jjs[j],(ind_diff-2):ind_diff]==0):
                                    new_diff_list.append(ind_diff)
                            DD['diff'] = new_diff_list
                        DD['ind'] = DD['diff'][np.argmin(np.abs(DD['diff']-DD['ind']))]
                    else:
                        DD['ind'] = DD['diff'][0]

            WD_ind = np.min([WD_rho['ind'],WD_u['ind'],WD_v['ind']])
        
            if len(X5m['diff'])==0:
                X5m['ind'] = WD_ind-2
            
            if (WD_ind-X5m['ind'])<2:
                X5m['ind'] = WD_ind-2
                
            dye_01[old_nt+t,j] = np.nanmean(dye_01_0[t,:,jjs[j],int(X5m['ind']):int(WD_ind)])
            dye_02[old_nt+t,j] = np.nanmean(dye_02_0[t,:,jjs[j],int(X5m['ind']):int(WD_ind)])
            u0[old_nt+t,j] = np.nanmean(u00[t,:,jjs[j],int(X5m['ind']):int(WD_ind)])
            v0[old_nt+t,j] = np.nanmean(v00[t,:,jjs[j],int(X5m['ind']):int(WD_ind)])

    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1


var_list = ['lon_rho','lat_rho','mask_rho','lonshore','latshore','h0',\
'x_rho','y_rho','xshore','yshore','rshore',\
'iis','jjs',\
'buoylon','buoylat','bi','bj',\
'Dwave','Hwave','Lwave','zeta',\
'dye_01','dye_02','u0','v0','ot']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/nearshore_variables_wavebuoy_5m_2017–2019.p'
pickle.dump(D,open(outfn,'wb'))
