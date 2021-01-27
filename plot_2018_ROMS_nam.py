#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean as cmo
import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))
    
# Build a code to plot BC's generated by code
#
# I want to compare it with the BC file generated for 2018
# with the LV3 files for 2018
    
plt.close('all')
dir0 = '/data0/ebrasseale/NADB2018/Input/'
LV4_grid = dir0+'GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
LV4_nam_2018 = dir0+'roms_nam_LV4_20171221_20180110.nc'
# LV4_nam_2018 = dir0+'roms_nam_LV4_20180331_20180420.nc'
# LV4_nam_2018 = dir0+'roms_nam_LV4_20180818_20180907.nc'
# LV4_nam_2018 = dir0+'roms_nam_LV4_20181216_20190101.nc'

dir0_nam = '/home/x1wu/SDTJRE_2018/mfiles/NAM_data/LV3'
LV3_nam_2018 = dir0_nam+ 'roms_nam_LV3_20171221_20180120.nc'
# LV3_nam_2018 = dir0_nam+ 'roms_nam_LV3_20180321_20180420.nc'
# LV3_nam_2018 = dir0_nam+ 'roms_nam_LV3_20180818_20180917.nc'
# LV3_nam_2018 = dir0_nam+ 'roms_nam_LV3_20181216_20190115.nc'

dgrd = nc.Dataset(LV4_grid)
lonr_LV4 = dgrd['lon_rho'][:]
latr_LV4 = dgrd['lat_rho'][:]
maskr_LV4 = dgrd['mask_rho'][:]

dlv4 = nc.Dataset(LV4_nam_2018)
uwind_LV4 = dlv4['Uwind'][:]
vwind_LV4 = dlv4['Vwind'][:]
pair_LV4 = dlv4['pair'][:]

dlv3 = nc.Dataset(LV3_nam_2018)
lonr_LV3 = dlv3['lon'][:]
latr_LV3 = dlv3['lat'][:]
uwind_LV3 = dlv3['Uwind'][:]
vwind_LV3 = dlv3['Vwind'][:]
pair_LV3 = dlv3['pair'][:]

# now match time indexes between the two grids
date0 = datetime(2018,1,1)

wt3 = dlv3['wind_time'][:]
wt3_list = []
for t in wt3:
    date = datetime(1999,1,1)+timedelta(days=t)
    wt3_list.append(date)
t3 = wt3_list.index(date0)

wt4 = dlv4['wind_time'][:]
wt4_list = []
for t in wt4:
    date = datetime(1999,1,1)+timedelta(days=t)
    wt4_list.append(date)
t4 = wt4_list.index(date0)

fig=plt.figure(figsize=(10,14))
gs = GridSpec(3,2)

#plot LV3 on left
# plot u wind
vmin = -5
vmax = 5
ax0 = fig.add_subplot(gs[0,0])
p=ax0.pcolormesh(lonr_LV3,latr_LV3,uwind_LV3[t3,:],cmap='BrBG',vmin=vmin,vmax=vmax)
cbaxes = inset_axes(ax0, width="4%", height="40%", loc=4,bbox_transform=ax0.transAxes,bbox_to_anchor=(-0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('wind (ms-1)',rotation=90,labelpad=10,fontweight='bold')

ax0.set_ylabel('latitude')
ax0.set_xlabel('longitude')
labeltext= 'LV3 \n'+date0.strftime("%m/%d/%Y") + '\nUwind'
ax0.text(0.1,0.9,labeltext,transform=ax0.transAxes,fontweight='bold',va='top')

# plot LV3 v wind
ax1 = fig.add_subplot(gs[1,0])
ax1.pcolormesh(lonr_LV3,latr_LV3,vwind_LV3[t3,:],cmap='BrBG',vmin=vmin,vmax=vmax)
ax1.set_ylabel('latitude')
ax1.set_xlabel('longitude')
labeltext= 'LV3 \n'+date0.strftime("%m/%d/%Y") + '\nVwind'
ax1.text(0.1,0.9,labeltext,transform=ax1.transAxes,fontweight='bold',va='top')

# plot LV3 air pressure
ax2 = fig.add_subplot(gs[2,0])
ax2.pcolormesh(lonr_LV3,latr_LV3,pair_LV3[t3,:],cmap='RdPu',vmin=1013,vmax=1021)
ax2.set_ylabel('latitude')
ax2.set_xlabel('longitude')
labeltext= 'LV3 \n'+date0.strftime("%m/%d/%Y") + '\nPair'
ax2.text(0.1,0.9,labeltext,transform=ax2.transAxes,fontweight='bold',va='top')

# LV4 uwind
ax3 = fig.add_subplot(gs[0,1])
ax3.pcolormesh(lonr_LV4,latr_LV4,uwind_LV4[t4,:],cmap='BrBG',vmin=vmin,vmax=vmax)
ax3.set_ylabel('latitude')
ax3.set_xlabel('longitude')
labeltext= 'LV4 \n'+date0.strftime("%m/%d/%Y") + '\nUwind'
ax3.text(0.1,0.9,labeltext,transform=ax3.transAxes,fontweight='bold',va='top')

# LV4 vwind
ax4 = fig.add_subplot(gs[1,1])
ax4.pcolormesh(lonr_LV4,latr_LV4,vwind_LV4[t4,:,:],cmap='BrBG',vmin=vmin,vmax=vmax)
ax4.set_ylabel('latitude')
ax4.set_xlabel('longitude')
labeltext= 'LV4 \n'+date0.strftime("%m/%d/%Y") + '\nVwind'
ax4.text(0.1,0.9,labeltext,transform=ax4.transAxes,fontweight='bold',va='top')

# LV4 air pressure
ax5 = fig.add_subplot(gs[2,1])
ax5.pcolormesh(lonr_LV4,latr_LV4,pair_LV4[t4,:],cmap='RdPu',vmin=1013,vmax=1021)
ax5.set_ylabel('latitude')
ax5.set_xlabel('longitude')
labeltext= 'LV4 \n'+date0.strftime("%m/%d/%Y") + '\nPair'
ax5.text(0.1,0.9,labeltext,transform=ax5.transAxes,fontweight='bold',va='top')

for ax in ax0,ax1,ax2,ax3,ax4,ax5:
    ax.axis([lonr_LV4.min(),lonr_LV4.max(),latr_LV4.min(),latr_LV4.max()])
    ax.contour(lonr_LV4,latr_LV4,maskr_LV4,levels=[0.5],colors='k')

dlv3.close()
dlv4.close()
dgrd.close

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2018_NAM_LV3_LV4_'+date.strftime("%Y.%m.%d")+'.png'
plt.savefig(out_fn)
