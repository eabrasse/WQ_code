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

# Build a code to plot BC's generated by code
#
# I want to compare it with the BC file generated for 2018
# with the LV3 files for 2018
    
plt.close('all')

LV4_grid = '/data0/ebrasseale/NADB2018/Input/GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
LV4_BC_2018 = '/data0/ebrasseale/NADB2018/Input/BC_LV4_20171117_20180615_Nz10_dye.nc'
# goes from Nov 17 2017 to June 15 2018
LV3_BC_2018 = '/home/x1wu/SDTJRE_2018/LV3_RUNFILES/Run2018/ocean_his_LV3_EPA20172018_00027.nc'
#goes from Jan 1 2018 to Jan 16 2018

# LV3_BC_2018 = '/home/x1wu/SDTJRE_2018/LV3_RUNFILES/Run2018/ocean_his_LV3_EPA20172018_00037.nc'
# #goes from Jun 1 2018 to Jun 15 2018
# LV3_BC_2018 = '/home/x1wu/SDTJRE_2018/LV3_RUNFILES/Run2018/ocean_his_LV3_EPA20172018_00038.nc'
# #goes from Jun 15 2018 to Jan 30 2018
# LV3_BC_2018 = '/home/x1wu/SDTJRE_2018/LV3_RUNFILES/Run2018/ocean_his_LV3_EPA20172018_00051.nc'
# #goes from Dec 15 2018 to Dec 31 2018

dgrd = nc.Dataset(LV4_grid)
lonr_lv4 = dgrd['lon_rho'][:]
latr_lv4 = dgrd['lat_rho'][:]

dlv4 = nc.Dataset(LV4_BC_2018)

dlv3 = nc.Dataset(LV3_BC_2018)
lonr_lv3 = dlv3['lon_rho'][:]
latr_lv3 = dlv3['lat_rho'][:]

var_name = 'temp'
var_west_LV4 = dlv4[var_name+'_west'][:]
var_south_LV4 = dlv4[var_name+'_south'][:]
var_lv3 = dlv3[var_name]

# now match time indexes between the two grids
t3 = 0
ott = dlv3['ocean_time'][:]
ot = ott[t3]
date = datetime(1999,1,1)+timedelta(seconds=ot)
ot_days = ot/(24*60*60)

vart = dlv4[var_name+'_time'][:]
t_diff = np.abs(ot_days-vart)
print('t_diff min = {:0.4f}'.format(t_diff.min()))
t4 = np.argmin(t_diff)

# extract relevant LV3 data along boundary
# unfortunately, grids are tilted and not plaid
# have to cycle through indexes in LV4 and find
# corresponding LV3 ones
#
# since LV4 grid is higher resolution,
# I might be able to reduce to "unique" LV3 indexes
ny4, nx4 = lonr_lv4.shape

ji_4to3_west = []

# build up western edge
for j in range(ny4):
    lon4 = lonr_lv4[j,0]
    lat4 = latr_lv4[j,0]
    # find nearest lv3 ind
    londiff = np.abs(lonr_lv3-lon4)
    latdiff = np.abs(latr_lv3-lat4)
    # grid of the radial difference from each grid point to the reference
    lonlatdiff = np.sqrt(londiff**2+latdiff*2)
    # take grid cell w/ smallest distance
    j_4to3_west = np.where(lonlatdiff==lonlatdiff.min())[0][0]
    i_4to3_west = np.where(lonlatdiff==lonlatdiff.min())[1][0]
    ji_4to3_west = ji_4to3_west + [[j_4to3_west, i_4to3_west]]

#reduce index pairs to only unique pairs on LV3 grid
LV3_ji_west = []
for ji in ji_4to3_west:
    if ji not in LV3_ji_west:
        LV3_ji_west = LV3_ji_west + [ji]

ny3 = len(LV3_ji_west)
nz3 = var_lv3.shape[1] # should have t,z,y,x indexes
var_lv3_west = np.zeros((nz3,ny3))
for ji in range(ny3):
    j,i = LV3_ji_west[ji]
    var_lv3_west[:,ji] = var_lv3[t3,:,j,i]

#repeat for southern edge
ji_4to3_south = []

for i in range(nx4):
    lon4 = lonr_lv4[0,i]
    lat4 = latr_lv4[0,i]
    # find nearest lv3 ind
    londiff = np.abs(lonr_lv3-lon4)
    latdiff = np.abs(latr_lv3-lat4)
    # grid of the radial difference from each grid point to the reference
    lonlatdiff = np.sqrt(londiff**2+latdiff*2)
    # take grid cell w/ smallest distance
    j_4to3_south = np.where(lonlatdiff==lonlatdiff.min())[0][0]
    i_4to3_south = np.where(lonlatdiff==lonlatdiff.min())[1][0]
    ji_4to3_south = ji_4to3_south + [[j_4to3_south, i_4to3_south]]

#reduce index pairs to only unique pairs on LV3 grid
LV3_ji_south = []
for ji in ji_4to3_south:
    if ji not in LV3_ji_south:
        LV3_ji_south = LV3_ji_south + [ji]

nx3 = len(LV3_ji_south)

var_lv3_south = np.zeros((nz3,nx3))
for ji in range(nx3):
    j,i = LV3_ji_south[ji]
    var_lv3_south[:,ji] = var_lv3[t3,:,j,i]


fig=plt.figure(figsize=(12,10))
gs = GridSpec(3,2)

#plot LV3 on left
# start with west boundary
ax0 = fig.add_subplot(gs[0,0])
p=ax0.pcolormesh(var_lv3_west,cmap='YlOrRd',vmin=5,vmax=20)
cbaxes = inset_axes(ax0, width="4%", height="40%", loc=4,bbox_transform=ax0.transAxes,bbox_to_anchor=(-0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('temp (C)',rotation=90,labelpad=10,fontweight='bold')

ax0.set_ylabel('vertical index')
ax0.set_xlabel('horizontal index')
labeltext= 'LV3 /n'+date.strftime("%m/%d/%Y") + '/n'+ var_name + '/n' + 'west'
ax0.text(0.1,0.9,labeltext,transform=ax0.transAxes,fontweight='bold')

# next LV3 south boundary
ax1 = fig.add_subplot(gs[1,0])
ax1.pcolormesh(var_lv3_south,cmap='YlOrRd',vmin=5,vmax=20)
ax1.set_ylabel('vertical index')
ax1.set_xlabel('horizontal index')
labeltext= 'LV3 /n'+date.strftime("%m/%d/%Y") + '/n'+ var_name + '/n' + 'south'
ax1.text(0.1,0.9,labeltext,transform=ax1.transAxes,fontweight='bold')


#finally add time series at bottom of western boundary
ax2 = fig.add_subplot(gs[2,:])
dt_list = []
nt3 = var_lv3.shape[0]
for tt in ott:
    dt_list.append(datetime(1999,1,1,0,0)+timedelta(seconds=tt))
j0,i0 = LV3_ji_west[0]
var_lv3_timeseries = var_lv3[:,-1,j0,i0]

ax2.plot(dt_list,var_lv3_timeseries,color='cornflowerblue',label='LV3')
ax2.plot(dt_list,var_west_LV4[t4:t4+nt3,-1,0],color='orange',label='LV4')
ax2.legend()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax2.set_xlabel('Time')

# next LV4 west boundary
ax3 = fig.add_subplot(gs[0,1])
ax3.pcolormesh(var_LV4_west[t4,:,:],cmap='YlOrRd',vmin=5,vmax=20)
ax3.set_ylabel('vertical index')
ax3.set_xlabel('horizontal index')
labeltext= 'LV4 /n'+date.strftime("%m/%d/%Y") + '/n'+ var_name + '/n' + 'west'
ax3.text(0.1,0.9,labeltext,transform=ax3.transAxes,fontweight='bold')

# next LV4 south boundary
ax4 = fig.add_subplot(gs[1,1])
ax4.pcolormesh(var_LV4_south[t4,:,:],cmap='YlOrRd',vmin=5,vmax=20)
ax4.set_ylabel('vertical index')
ax4.set_xlabel('horizontal index')
labeltext= 'LV4 /n'+date.strftime("%m/%d/%Y") + '/n'+ var_name + '/n' + 'south'
ax4.text(0.1,0.9,labeltext,transform=ax4.transAxes,fontweight='bold')

dlv3.close()
dlv4.close()
dgrd.close

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2018_BC_'+var_name+'_LV3_LV4.png'
plt.savefig(out_fn)
