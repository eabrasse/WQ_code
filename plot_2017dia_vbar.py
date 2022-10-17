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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean as cmo

plt.close('all')
c10 = plt.get_cmap('tab10')
# useful plotting tools, colormap and test box props
landcol = 'lightgray'
seacol = 'white'

cmap_mask = matplotlib.colors.ListedColormap([landcol,seacol])

# home = '/dataSIO/ebrasseale/'

#build up list of file names for model output
dir0 = '/dataSIO/ebrasseale/NADB2017_0_dia/'

fn = dir0+'ocean_dia_NADB_0_new_00001.nc'

ds = nc.Dataset(fn)


var_list = ['vbar_cor','vbar_hadv','vbar_hjvf','vbar_wbrk','vbar_prsgrd','vbar_sstr','vbar_bstr','vbar_accel','vbar_hvisc']

# identify a few random places to plot
lonv = ds['lon_v'][:]
latv = ds['lat_v'][:]
maskv = ds['mask_v'][:]

# # z,y,x
# lat0,lon0 = 32.56957,-117.16880 # CDIP buoy location
# i0,j0 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon0,lat0)
# loc0 = [5,j0,i0]
#
# lat1,lon1 = 32.66,-117.13 # middle of SD Bay
# i1,j1 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon1,lat1)
# loc1 = [9,j1,i1]
#
# lat2,lon2 = 32.5,-117.2 # further south and offshore
# i2,j2 = wqfun.find_ll_inds(lon_rho,lat_rho,mask_rho,lon2,lat2)
# loc2 = [2,j2,i2]

# note: z's were chosen randomly

nvars = len(var_list)
nrows =2
if np.mod(nvars,nrows)>0:
    ncols = int(np.floor(nvars/nrows))+int(np.mod(nvars,nrows))
else:
    ncols = int(nvars/nrows)
fig,axs=plt.subplots(nrows=nrows,ncols=ncols,figsize=(12,8))

vmaxs = []
t=-1
for varname,ax in zip(var_list,axs.ravel()):
    
    # draw map
    ax.pcolormesh(lonv,latv,maskv,cmap=cmap_mask,shading='nearest',zorder=0)
    wqfun.dar(ax)
    
    # load in and mask data
    var = ds[varname][t,:]
    var = np.ma.masked_where(var==0,var)
    
    # label
    ax.text(0.9,0.9,varname,color='k',fontweight='bold',transform=ax.transAxes,ha='right',va='top')
    
    # try to determine what the best colorbar axis is going to be
    # it needs to be symmetric around zero
    vmaxs.append(np.max(np.abs(var)))

vmax = 0.8*np.max(vmaxs)

for ax in axs.ravel():
    # plot data
    p = ax.pcolormesh(lonv,latv,var,cmap=cmo.cm.balance,vmin=-vmax,vmax=vmax)

cbaxes = inset_axes(axs[-1,ncols-1], width="6%", height="40%", loc='lower left',bbox_transform=axs[-1,ncols-1].transAxes,bbox_to_anchor=(0.,0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.set_label('m/s2',fontsize=8)

for ax in axs[0,:]:
    ax.set_xlabel('Longitude')
for ax in axs[1:,:]:
    ax.set_xticklabels([''])
for ax in axs[:,0]:
    ax.set_ylabel('Latitude')
for ax in axs[:,1:].ravel():
    ax.set_yticklabels([''])
    
fig.subplots_adjust(left=0.08,right=0.98,top=0.95,wspace=0.3)
outfn = home+ f'WQ_plots/NADB2017_vbar_dia_snapshot.jpg'
plt.savefig(outfn)
plt.close()

# clean up
ds.close()