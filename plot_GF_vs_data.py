#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experime.
"""

# setup
import os
import sys
import wqfun
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import matplotlib.dates as mdates

plt.close('all')
c10 = plt.get_cmap('tab10')

home = '/dataSIO/ebrasseale/'

#build up list of file names for model output
dir0 = '/dataSIO/NADB2017/NADB2017_0_NEW/'
dir1 = '/dataSIO/NADB2018/'
dir2 = '/dataSIO/NADB2019/'

f_list = []
for my_dir in [dir0,dir1,dir2]:
    f_list0 = os.listdir(my_dir)
    f_list0.sort()
    f_list0 = [my_dir+x for x in f_list0 if x[:9]=='ocean_his']
    f_list.extend(f_list0)
    
#choose files to plot
f_list_inds = [0, 50, -1]
f_list1 = [f_list[ind] for ind in f_list_inds]
nfile = len(f_list1)

ds_list = [nc.Dataset(fn) for fn in f_list1]

# find start and end times of each dataset
t0_list = [ds['ocean_time'][:][0] for ds in ds_list]
t1_list = [ds['ocean_time'][:][-1] for ds in ds_list]

# make list of datetime lists: len(dt_list) = nfile, each item in dt_list is a list of ~240 datetimes
dt_list = [[datetime(1999,1,1)+timedelta(seconds=t) for t in ds['ocean_time'][:][:]] for ds in ds_list]

# load in Godin Filtered data set
GF_fn = home+'WQ_data/ocean_daily_gf_NADB2017-2018-2019.nc'
dsgf = nc.Dataset(GF_fn)

# make datetime list of GF dataset
dt_listGF = [datetime(1999,1,1)+timedelta(seconds=t) for t in dsgf['ocean_time'][:]]

# find start and end indeces for the godin-filtered dataset for each file
gft0 = [np.argwhere((dsgf['ocean_time'][:]-t0)>0)[0][0] for t0 in t0_list]
gft1 = [np.argwhere((dsgf['ocean_time'][:]-t1)<0)[-1][0] for t1 in t1_list]

testing=True
if testing:
    # just plot one variable
    var_time_list = ['zeta']
else:
    # make a list of time-varying variables
    var_time_list = [v_name for v_name,varin in dsgf.variables.items() if 'ocean_time' in varin.dimensions]

# identify a few random places to plot
lon_rho = dsgf['lon_rho'][:]
lat_rho = dsgf['lat_rho'][:]
mask_rho = dsgf['mask_rho'][:]

# z,y,x
lat0,lon0 = 32.56957,-117.16880 # CDIP buoy location
i0,j0 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon0,lat0)
loc0 = [5,j0,i0]

lat1,lon1 = 32.66,-117.13 # middle of SD Bay
i1,j1 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon1,lat1)
loc1 = [9,j1,i1]

lat2,lon2 = 32.5,-117.2 # further south and offshore
i2,j2 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon2,lat2)
loc2 = [2,j2,i2]

# note: z's were chosen randomly

loc_list = [loc0,loc1,loc2]
nloc = len(loc_list)

for vname in var_time_list:
    
    # no need to plot ocean time
    # all other variables have at least 3 dimensions (time+x+y)
    ndim = len(dsgf[vname].dimensions)
    if ndim<3:
        continue
    
    # open figure for this variable
    fig=plt.figure(figsize=(10,6))
    gs = GridSpec(nloc,nfile+1)
    
    # draw map on left for extraction location reference
    ax_map = fig.add_subplot(gs[:,0])
    ax_map.contour(lon_rho,lat_rho,mask_rho,levels=[0.5],colors='gray',label=None)
    count = 0
    for loc in loc_list:
        ax_map.plot(lon_rho[loc[1],loc[2]],lat_rho[loc[1],loc[2]],marker='*',color='none',markersize=15,mec='k',mfc=c10(count))
        count+=1
    wqfun.dar(ax_map)
    ax_map.set_xlabel('longitude')
    ax_map.set_ylabel('latitude')
    ax_map.set_title('Extraction locations')
    
    # load in godin-filtered variable
    var_GF = dsgf[vname][:]
    
    tcount =0
    for ds in ds_list:
        var = ds[vname][:]
        
        lcount=0
        for loc in loc_list:
            
            # identify indexes for this location
            if ndim==3:
                # var dimensions are t, y, x
                inds_list = loc[1:] # loc has [z, y, x], just take [y, x]
                inds_list.insert(0,slice(None)) # include all t
                inds = tuple(inds_list) # convert list to tuple for indexing
                
                gf_inds_list = loc[1:] # loc has [z, y, x], just take [y, x]
                gf_inds_list.insert(0,slice(gft0[tcount],gft1[tcount])) # t index from gft0 to gft1
                gf_inds = tuple(gf_inds_list) # convert list to tuple for indexing
                
            elif ndim==4:
                # var dimensions are t, z, y, x
                inds_list = loc.copy() # loc has [z, y, x], take all
                inds_list.insert(0,slice(None)) # include all t
                inds = tuple(inds_list) # convert list to tuple for indexing
                
                gf_inds_list = loc.copy # loc has [z, y, x], take all
                gf_inds_list.insert(0,slice(gft0[tcount],gft1[tcount])) # t index from gft0 to gft1
                gf_inds = tuple(gf_inds_list) # convert list to tuple for indexing
            
            # generate time series just at local point
            var_local = var[inds]
            var_GF_local = var_GF[gf_inds]
            
            # godin filter local time series to test
            var_hourly_gf_local = wqfun.filt_godin(var_local)
            
            # grab subplot handle
            ax = fig.add_subplot(gs[lcount,tcount+1])
            
            # plot unfilted local data
            ax.plot(dt_list[tcount],var_local,lw=1.5,ls='solid',color=c10(lcount),label='data')
            
            # plot locally filtered data
            ax.plot(dt_list[tcount],var_hourly_gf_local,lw=1.0,ls='dotted',color=c10(lcount),label='1D GF data')
            
            # plot globally filtered and subsampled local data
            ax.plot(dt_listGF[gft0[tcount]:gft1[tcount]],var_GF_local,marker='x',ls='None',color='k',markersize=8,label='multi-D GF, 24hr subsampled data')
            
            if (tcount==0) and (lcount==0):
                ax.legend()
                
            # add units to y-axis on far left hand side only
            if tcount<1:
                ax.set_ylabel(ds[vname].units)
            
            # add time to x-axis on bottom plots only
            if lcount<nloc-1:
                ax.set_xticklabels([''])
            else:
                
                ax.set_xlabel('date')
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d, %Y"))
                plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
            
            lcount+=1
        tcount+=1
    
    fig.subplots_adjust(left=0.05,right=0.95,top=0.85)
    outfn = home+ f'WQ_plots/ocean_daily_gf_NADB2017-2018-2019_compare_{vname}.jpg'
    plt.savefig(outfn)
    plt.close()

# clean up
dsgf.close()
for ds in ds_list:
    ds.close()