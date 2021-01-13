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
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# code to compare Punta Bandera river input file with model behavior
# Qin should = Qout, dye_in should = dye_out, etc.

home = '/data0/ebrasseale/'

dir0 = '/data0/NADB2017/NADB2017_0/Output/'

f_list = os.listdir(dir0)
f_list.sort()
riv_fn = [x for x in f_list if x[0]=='r'][0]
f_list = [x for x in f_list if x[:17]=='ocean_his_NADB_0_']

testing=True
if testing:
    f_list = f_list[0]

# nfiles = len(f_list)
dsr = nc.Dataset(dir0+riv_fn)
i00 = dsr['river_Xposition'][5]-1
j00 = dsr['river_Eposition'][5]

i00 = int(i00)
j00 = int(j00)

i0 = int(i00-5)
j0 = int(j00-5)

i1 = int(i00+5+1)
j1 = int(j00+5+1)

# count=0
# for fn in f_list:
ds = nc.Dataset(dir0+f_list)
# if count == 0:
lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
lon_u = ds['lon_u'][:]
lat_u = ds['lat_u'][:]
mask_rho = ds['mask_rho'][:]
masked_rho = 1 - mask_rho
wetdry_mask_rho = ds['wetdry_mask_rho'][:]
wetdry_masked_rho = 1 - wetdry_mask_rho

#variables to plot
dye_01 = ds['dye_01'][:]
salt = ds['salt'][:]
nt,nz,ny,nx = salt.shape

for t in range(2): # change to range(nt) when ready to go
    # initialize figure and axis
    fig = plt.figure(figsize=(14,6))
    axdye = fig.add_subplot(1,2,1)
    
    # plot dye
    pd=axdye.pcolormesh(lon_rho[j0:j1,i0:i1],lat_rho[j0:j1,i0:i1],dye_01[t,-1,j0:j1,i0:i1],cmap='Greens',linewidths=1,shading='nearest')
    
    # add colorbar
    cbaxes = inset_axes(axdye, width="4%", height="80%", loc=4,bbox_transform=axdye.transAxes,bbox_to_anchor=(0.075,0.,1,1))
    cb = fig.colorbar(pd, cax=cbaxes, orientation='vertical')
    cb.set_label('dye conc')
    
    # add other labels
    axdye.set_title('Surface dye')
    axdye.set_xlabel('Longitude')
    axdye.set_ylabel('Latitude')
    
    axdye.text(0.1,0.2,'x = land mask',transform=axdye.transAxes,color='black',va='center',ha='left',fontweight='bold')
    axdye.text(0.1,0.1,'o = wetdry mask',transform=axdye.transAxes,color='magenta',va='center',ha='left',fontweight='bold')
    
    # initialize next axis
    axsalt = fig.add_subplot(1,2,2)
    
    #plot salt
    ps=axsalt.pcolormesh(lon_rho[j0:j1,i0:i1],lat_rho[j0:j1,i0:i1],salt[t,-1,j0:j1,i0:i1],cmap='Blues',linewidths=1,shading='nearest')
    
    # add colorbar
    cbaxes = inset_axes(axsalt, width="4%", height="80%", loc=4,bbox_transform=axsalt.transAxes,bbox_to_anchor=(0.075,0.,1,1))
    cb = fig.colorbar(ps, cax=cbaxes, orientation='vertical')
    cb.set_label('salinity')
    
    # add other labels
    axsalt.set_title('Surface salinity')
    axsalt.set_xlabel('Longitude')
    axsalt.set_ylabel('Latitude')
    
    # add context details to both plots
    for ax in axdye, axsalt:
        # add river indicator
        ax.plot(lon_u[j00,i00],lat_u[j00,i00],marker='<',markeredgecolor='black',markerfacecolor='yellow',markersize=8)
        
        # scatter land mask
        ax.scatter(lon_rho[j0:j1,i0:i1],lat_rho[j0:j1,i0:i1],10*masked_rho[j0:j1,i0:i1],c='None',marker='x',edgecolors='black')
        
        # contour wet-dry mask
        ax.scatter(lon_rho[j0:j1,i0:i1],lat_rho[j0:j1,i0:i1],10*wetdry_masked_rho[t,j0:j1,i0:i1],c='None',marker='o',edgecolors='magenta')
    
    #save, close, start next
    # plt.tight_layout()
    outfn = home+'WQ_plots/river_zoom_movie/fig_{:04}.png'.format(t)
    plt.savefig(outfn)
    plt.close()
