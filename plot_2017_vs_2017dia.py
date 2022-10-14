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
dir1 = '/dataSIO/ebrasseale/NADB2017_0_dia/'

fn_xw = dir0+'ocean_his_NADB_0_new_00001.nc'
fn_dia = dir1+'ocean_his_NADB_0_new_00001.nc'


ds_xw = nc.Dataset(fn_xw)
ds_dia = nc.Dataset(fn_dia)

testing=True
if testing:
    # just plot one variable
    var_time_list = ['zeta','u','v','w','salt','temp']
else:
    # make a list of time-varying variables
    var_time_list = [v_name for v_name,varin in dsgf.variables.items() if 'ocean_time' in varin.dimensions]

# identify a few random places to plot
lon_rho = ds_xw['lon_rho'][:]
lat_rho = ds_xw['lat_rho'][:]
mask_rho = ds_xw['mask_rho'][:]

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
    ndim = len(ds_xw[vname].dimensions)
    if ndim<3:
        continue
    
    # open figure for this variable
    fig=plt.figure(figsize=(12,8))
    gs = GridSpec(nloc,2)
    
    # draw map on left for extraction location reference
    ax_map = fig.add_subplot(gs[:,0])
    ax_map.contour(lon_rho,lat_rho,mask_rho,levels=[0.5],colors='gray')
    count = 0
    for loc in loc_list:
        ax_map.plot(lon_rho[loc[1],loc[2]],lat_rho[loc[1],loc[2]],marker='*',color='none',markersize=15,mec='k',mfc=c10(count))
        count+=1
    wqfun.dar(ax_map)
    ax_map.set_xlabel('longitude')
    ax_map.set_ylabel('latitude')
    ax_map.set_title('Extraction locations')
        
    lcount=0
    for loc in loc_list:
        
        # identify indexes for this location
        if ndim==3:
            # var dimensions are t, y, x
            inds_list = loc[1:] # loc has [z, y, x], just take [y, x]
            inds_list.insert(0,slice(None)) # include all t
            inds = tuple(inds_list) # convert list to tuple for indexing
            
        elif ndim==4:
            # var dimensions are t, z, y, x
            inds_list = loc.copy() # loc has [z, y, x], take all
            inds_list.insert(0,slice(None)) # include all t
            inds = tuple(inds_list) # convert list to tuple for indexing
        
        
        # grab subplot handle
        ax = fig.add_subplot(gs[lcount,1])
        
        # plot unfilted local data
        var_xw = ds_xw[vname][:]
        ax.plot(dt_list,var_xw[inds],lw=1.0,ls='solid',color=c10(0))
        var_dia = ds_dia[vname][:]
        ax.plot(dt_list,var_dia[inds],lw=1.0,ls='dashed',color=c10(1))
        
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
        
        # add units to y-axis on far left hand side only
        try:
            ax.set_ylabel(ds[vname].units)
        except:
            ax.set_ylabel(vname)
        
        # add time to x-axis on bottom plots only
        if lcount<nloc-1:
            ax.set_xticklabels([''])
        else:
            
            ax.set_xlabel('date')
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %-d, %Y"))
            plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
        
        lcount+=1
    
    fig.subplots_adjust(left=0.08,right=0.98,top=0.95,wspace=0.3)
    outfn = home+ f'WQ_plots/NADB2017_dia_validation/compare_NADB2017_NADB2017dia_{vname}.jpg'
    plt.savefig(outfn)
    plt.close()

# clean up
ds_xw.close()
ds_dia.close()