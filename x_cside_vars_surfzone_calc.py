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

# this piece of code is to test extraction of variables from the COAWST model
# using a variable surf zone width
# it only extracts waves at one location (in line with NDBC buoy Station 46235 - Imperial Beach Nearshore, CA (155))
# to simulate how I expect the model to be adapted for other users

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


testing=True
if testing:
    f_list = f_list[:1]

# ref_depth = 5
buoylat = 32.570
buoylon = -117.169

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
        mask_diff = {}
        lonshore = np.zeros((len(jjs)))
        latshore = np.zeros((len(jjs)))
        for j in range(len(jjs)):
            mask_diff[j] = np.where(np.diff(mask_rho[jjs[j],:]))[0]
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

# derived values
hb = np.zeros((NT)) # depth of wave breaking
Lsz = np.zeros((NT,nj)) # width of surf zone

# time
ot = np.zeros((NT)) # time

# Now do the extraction and processing
tt=0
old_nt = 0
for fn in f_list:
    print('file {:d} of {:d}'.format(tt,nfiles))
    ds = nc.Dataset(fn)

    # select wave direction and significant wave height
    wetdry_mask_rho = ds['wetdry_mask_rho'][:]
    wetdry_mask_u = ds['wetdry_mask_u'][:]
    wetdry_mask_v = ds['wetdry_mask_v'][:]
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
    
    # calculate wave breaking depth based on wave height at buoy
    hb0 = Hwave0[:,bj,bi]/0.78 # from Hwave = gamma * depth, gamma = 0.78 from McCowan (1984) [section 5.7.3 of Kumar CEE 473
    # because the land value of H is 0.25 m...
    hb0[hb0<0.25] = 0.26
    hb[old_nt:old_nt+nt] = hb0

    for j in range(nj):
        #loop through time steps, 
        # because extraction indexes depend on time-varying wetdry mask
        for t in range(nt):
            # find the edge of the mask
            wd_mask_diff_rho = np.where(np.diff(wetdry_mask_rho[t,jjs[j],:]))[0]
            wd_mask_diff_u = np.where(np.diff(wetdry_mask_u[t,jjs[j],:]))[0]
            wd_mask_diff_v = np.where(np.diff(wetdry_mask_v[t,jjs[j],:]))[0]
            
            
            #find where depth crosses from deeper than ref_depth to shallower
            depth_diff = np.where(np.diff(np.sign(H[t,jjs[j],:]-hb0[t])))[0]
    
            #if multiple edges, north of TJRE
            if (len(mask_diff[j])>1)&(lat_rho[jjs[j],0]>32.6):
                #look for the edge closest to the previously identified edge
                x_wd_ind_rho = wd_mask_diff_rho[np.argmin(np.abs(x_wd_ind_rho-wd_mask_diff_rho))]
                x_wd_ind_u = wd_mask_diff_u[np.argmin(np.abs(x_wd_ind_u-wd_mask_diff_u))]
                x_wd_ind_v = wd_mask_diff_v[np.argmin(np.abs(x_wd_ind_v-wd_mask_diff_v))]
                x_sz_ind = depth_diff[np.argmin(np.abs(x_5m_ind-depth_diff))]

            #if multiple edges, south of TJRE
            elif (len(mask_diff[j])>1)&(lat_rho[jjs[j],0]<32.6):
                #do outermost edge
                x_wd_ind_rho = wd_mask_diff_rho[0]
                x_wd_ind_u = wd_mask_diff_u[0]
                x_wd_ind_v = wd_mask_diff_v[0]
                x_sz_ind = depth_diff[0]

            elif len(mask_diff[j])==1:
                x_wd_ind_rho = wd_mask_diff_rho[0]
                x_wd_ind_u = wd_mask_diff_u[0]
                x_wd_ind_v = wd_mask_diff_v[0]
                x_sz_ind = depth_diff[0]
                
            x_wd_ind = np.min([x_wd_ind_rho,x_wd_ind_u,x_wd_ind_v])

            if (x_wd_ind-x_sz_ind)<2:
                x_sz_ind = x_wd_ind-2
                
            
            
            #go offshore of the wet/dry mask by a tad
            # x_wd_ind = x_wd_ind - 2
            dye_01[old_nt+t,j] = np.nanmean(dye_01_0[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
            dye_02[old_nt+t,j] = np.nanmean(dye_02_0[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
            u0[old_nt+t,j] = np.nanmean(u00[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
            v0[old_nt+t,j] = np.nanmean(v00[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
            
            # calculate surfzone width
            Lsz_x = x_rho[jjs[j],int(x_wd_ind)]-x_rho[jjs[j],int(x_sz_ind)]
            Lsz_y = y_rho[jjs[j],int(x_wd_ind)]-y_rho[jjs[j],int(x_sz_ind)]
            Lsz[old_nt+t,j] = np.sqrt(Lsz_x**2+Lsz_y**2)
    
    ot[old_nt:old_nt+nt] = ocean_time
    old_nt += nt
    
    ds.close()
    tt+=1


var_list = ['lon_rho','lat_rho','mask_rho','lonshore','latshore','h0',\
'x_rho','y_rho','xshore','yshore','rshore',\
'iis','jjs','shoreangle',\
'buoylon','buoylat','bi','bj'\
'Dwave','Hwave','Lwave','zeta',\
'Lsz','hb','ot'\
'dye_01','dye_02','u0','v0']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/shoreline_variables_SZ_2017–2018.p'
pickle.dump(D,open(outfn,'wb'))
