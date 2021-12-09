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
import wqfun


home = '/data0/ebrasseale/'

ind_fn = home+'WQ_data/shore_buoy_inds.p'
Dinds = pickle.load(open(ind_fn,'rb'))
for var in Dinds.keys():
    locals()[var]=Dinds[var]

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

fn = f_list[0]
ds = nc.Dataset(dir0+fn)

lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
mask_rho = ds['mask_rho'][:]

ot = ds['ocean_time'][:]
NT = np.shape(ot)
nj = len(jjs)

lonshore = np.zeros((nj))
latshore = np.zeros((nj))
for j in range(nj):
    mask_diff[j] = np.where(np.diff(mask_rho[jjs[j],:]))[0]
    lonshore[j] = lon_rho[jjs[j],iis[j]]
    latshore[j] = lat_rho[jjs[j],iis[j]]

x_rho,y_rho = wqfun.ll2xy(lon_rho,lat_rho,lon_rho.min(),lat_rho.min())


buoylat = 32.570
buoylon = -117.169
bi, bj = wqfun.find_nearest_ind_2D(lon_rho,lat_rho,buoylon,buoylat)

# u0 = np.zeros((NT,nj)) # velocity in one direction (eight degrees from E-W)
# v0 = np.zeros((NT,nj)) # velocity in an orthogonal direction (eight degreees from N-S)

# derived values
hb = np.zeros((NT)) # depth of wave breaking
Lsz = np.zeros((NT,nj)) # width of surf zone

lon_sz = np.zeros((NT,nj))
lat_sz = np.zeros((NT,nj))
lon_wd = np.zeros((NT,nj))
lat_wd = np.zeros((NT,nj))


#imitating the way it's done in the code
wetdry_mask_rho = ds['wetdry_mask_rho'][:]
wetdry_mask_u = ds['wetdry_mask_u'][:]
wetdry_mask_v = ds['wetdry_mask_v'][:]

# Dwave0 = ds['Dwave'][:]
Hwave0 = ds['Hwave'][:]
# Lwave0 = ds['Lwave'][:]
zeta0 = ds['zeta'][:]

u00 = ds['u'][:]
v00 = ds['v'][:]

H = h0+zeta0

#the wave extractions are from the same points in space no matter what the time
#so they can exist outside of the shoreline indexing loop
# Dwave = Dwave0[:,bj,bi]
Hwave = Hwave0[:,bj,bi]
# Lwave = Lwave0[:,bj,bi]
zeta = zeta0[:,bj,bi]

# calculate wave breaking depth based on wave height at buoy
hb = Hwave0[:,bj,bi]/0.78 # from Hwave = gamma * depth, gamma = 0.78 from McCowan (1984) [section 5.7.3 of Kumar CEE 473
# because the land value of H is 0.25 m...
hb[hb<0.25] = 0.26

count_wd_mask_diff_rho = 0
count_wd_mask_diff_u = 0
count_wd_mask_diff_v = 0

count_depth_diff = 0

for t in range(NT):
    #loop through time steps, 
    # because extraction indexes depend on time-varying wetdry mask
    for j in range(nj):
        # find the edge of the mask
        wd_mask_diff_rho = np.where(np.diff(wetdry_mask_rho[t,jjs[j],:iis[j]]))[0]
        wd_mask_diff_u = np.where(np.diff(wetdry_mask_u[t,jjs[j],:iis[j]]))[0]
        wd_mask_diff_v = np.where(np.diff(wetdry_mask_v[t,jjs[j],:iis[j]]))[0]
        
        
        #find where depth crosses from deeper than ref_depth to shallower
        depth_diff = np.where(np.diff(np.sign(H[t,jjs[j],:iis[j]]-hb[t])))[0]
        
        if len(wd_mask_diff_rho)>1:
            count_wd_mask_diff_rho += 1
        if len(wd_mask_diff_u)>1:
            count_wd_mask_diff_u += 1
        if len(wd_mask_diff_v)>1:
            count_wd_mask_diff_v += 1
        if len(depth_diff)>1:
            count_depth_diff += 1

        # #if multiple edges, north of TJRE
        # if (len(mask_diff[j])>1)&(lat_rho[jjs[j],0]>32.6):
        #     #look for the edge closest to the previously identified edge
        #     x_wd_ind_rho = wd_mask_diff_rho[np.argmin(np.abs(x_wd_ind-wd_mask_diff_rho))]
        #     x_wd_ind_u = wd_mask_diff_u[np.argmin(np.abs(x_wd_ind-wd_mask_diff_u))]
        #     x_wd_ind_v = wd_mask_diff_v[np.argmin(np.abs(x_wd_ind-wd_mask_diff_v))]
        #     x_sz_ind = depth_diff[np.argmin(np.abs(x_sz_ind-depth_diff))]
        #
        # #if multiple edges, south of TJRE
        # elif (len(mask_diff[j])>1)&(lat_rho[jjs[j],0]<32.6):
        #     #do outermost edge
        #     x_wd_ind_rho = wd_mask_diff_rho[0]
        #     x_wd_ind_u = wd_mask_diff_u[0]
        #     x_wd_ind_v = wd_mask_diff_v[0]
        #     x_sz_ind = depth_diff[0]
        #
        # elif len(mask_diff[j])==1:
        
        x_wd_ind_rho = wd_mask_diff_rho[0]
        x_wd_ind_u = wd_mask_diff_u[0]
        x_wd_ind_v = wd_mask_diff_v[0]
        x_sz_ind = depth_diff[0]
            
        x_wd_ind = np.min([x_wd_ind_rho,x_wd_ind_u,x_wd_ind_v])

        if (x_wd_ind-x_sz_ind)<2:
            x_sz_ind = x_wd_ind-2
            

        # u0[t,j] = np.nanmean(u00[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
        # v0[t,j] = np.nanmean(v00[t,:,jjs[j],int(x_sz_ind):int(x_wd_ind)])
        
        # calculate surfzone width
        Lsz_x = x_rho[jjs[j],int(x_wd_ind)]-x_rho[jjs[j],int(x_sz_ind)]
        Lsz_y = y_rho[jjs[j],int(x_wd_ind)]-y_rho[jjs[j],int(x_sz_ind)]
        Lsz[t,j] = np.sqrt(Lsz_x**2+Lsz_y**2)
        
        lon_sz[t,j] = lon_rho[jjs[j],iis[j]]
        lat_sz[t,j] = lat_rho[jjs[j],iis[j]]
        lon_wd[t,j] = lon_rho[jjs[j],iis[j]]
        lat_wd[t,j] = lat_rho[jjs[j],iis[j]]


    fig = plt.figure(figsize=(6,8))
    ax = fig.gca()
    ax.contour(lon_rho,lat_rho,mask_rho,colors='k',levels=[1],linewidths=0.5,alpha=1.0)
    ax.scatter(lon_sz[t,:],lat_sz[t,:],marker='>',c='m')
    ax.scatter(lon_wd[t,:],lat_wd[t,:],marker='<',c='g')
    ax.set_ylim([latshore.min(),latshore.max()])
    ax.set_xlim([[-117.225,-117.075]])
    wqfun.dar(ax)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')

    outfn = home + f'WQ_plots/shoreline_wd_sz/figure{t:03}.png'
    plt.savefig(outfn)
    plt.close('all')
ds.close()

var_list = ['lon_rho','lat_rho','mask_rho','lonshore','latshore','h0',\
'lon_sz','lat_sz','lon_wd','lat_wd',\
'iis','jjs','shoreangle',\
'buoylon','buoylat','bi','bj',\
'Hwave','zeta',\
'Lsz','hb','ot']

D = dict()
for var in var_list:
    D[var]=locals()[var]

outfn = home + 'WQ_data/wetdry_surfzone_locations.p'
pickle.dump(D,open(outfn,'wb'))

