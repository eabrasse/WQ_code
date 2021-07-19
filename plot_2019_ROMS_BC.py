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

# LV4_grid = '/data0/NADB2019/Input/GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
LV4_grid = '/data0/NADB2018/ocean_his_NADB2018_00001.nc'
LV4_BC_2019 = '/data0/NADB2019/Input/BC_LV4_20181227_20190625_Nz10_dye.nc'
# goes from Dec 27 2018 to June 25 2019
# LV4_BC_2019= '/data0/NADB2019/Input/BC_LV4_20190610_20200106_Nz10_dye.nc'
# goes from Jun 10 2019 to Jan 6 2020

# LV3_BC_2019 = '/data0/SDTJRE_LV3/LV3_2016_2020/ocean_his_LV3_EXT_20172018_00051.nc'
#goes from Dec 27 2018 to Jan 11 2019
LV3_BC_2019 = '/data0/SDTJRE_LV3/LV3_2016_2020/ocean_his_LV3_EXT_20172018_00062.nc'
# #goes from Jun 10 2019 to Jun 25 2019
# LV3_BC_2019 = '/data0/SDTJRE_LV3/LV3_2016_2020/ocean_his_LV3_EXT_20172018_00063.nc'
# #goes from Jun 25 2019 to Jul 10 2019
# LV3_BC_2019 = '/data0/SDTJRE_LV3/LV3_2016_2020/ocean_his_LV3_EXT_20172018_00075.nc'
# #goes from Dec 22 2019 to Jan 6 2020

dgrd = nc.Dataset(LV4_grid)
lonr_lv4 = dgrd['lon_rho'][:]
latr_lv4 = dgrd['lat_rho'][:]
sr_lv4 = dgrd['s_rho'][:]

dlv4 = nc.Dataset(LV4_BC_2019)

dlv3 = nc.Dataset(LV3_BC_2019)
lonr_lv3 = dlv3['lon_rho'][:]
latr_lv3 = dlv3['lat_rho'][:]
sr_lv3 = dlv3['s_rho'][:]
maskr_lv3 = dlv3['mask_rho'][:]

Dsalt = {'var_name':'salt','axlabel':'salinity (psu)','vmin':33.6,'vmax':33.8}
Dtemp = {'var_name':'temp','axlabel':'temp (C)','vmin':16.5,'vmax':18}

Dvar = Dtemp
var_name = Dvar['var_name'][:]
var_west_LV4 = dlv4[var_name+'_west'][:]
var_south_LV4 = dlv4[var_name+'_south'][:]
var_lv3 = dlv3[var_name]

# now match time indexes between the two grids
t30 = 0
ott = dlv3['ocean_time'][:]
ot0 = ott[t30]
date = datetime(1999,1,1)+timedelta(seconds=ot0)
ot0_days = ot0/(24*60*60)

vart = dlv4[var_name+'_time'][:]
t_diff = np.abs(ot0_days-vart)
print('t_diff min = {:0.4f}'.format(t_diff.min()))
t40 = np.argmin(t_diff)

# t3 = -1
# ot = ott[t3]
# date = datetime(1999,1,1)+timedelta(seconds=ot)
# ot_days = ot/(24*60*60)
#
# t_diff = np.abs(ot_days-vart)
# t4 = np.argmin(t_diff)
t3 = t30
t4 = t40

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
unique_j = []
LV3_ji_west = []
for ji in ji_4to3_west:
    j,i = ji
    # if ji not in LV3_ji_west:
    if j not in unique_j:
        unique_j = unique_j + [j]
        LV3_ji_west = LV3_ji_west + [ji]

ny3 = len(LV3_ji_west)
nz3 = var_lv3.shape[1] # should have t,z,y,x indexes
nz4 = var_west_LV4.shape[1]
var_lv3_west = np.zeros((nz3,ny3))
lat_lv3_west = np.zeros((ny3))
for ji in range(ny3):
    j,i = LV3_ji_west[ji]
    var_lv3_west[:,ji] = var_lv3[t3,:,j,i]
    lat_lv3_west[ji] = latr_lv3[j,i]

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
unique_i = []
LV3_ji_south = []
for ji in ji_4to3_south:
    j,i = ji
    # if ji not in LV3_ji_south:
    if i not in unique_i:
        unique_i = unique_i + [i]
        LV3_ji_south = LV3_ji_south + [ji]

nx3 = len(LV3_ji_south)

var_lv3_south = np.zeros((nz3,nx3))
lon_lv3_south = np.zeros((nx3))
for ji in range(nx3):
    j,i = LV3_ji_south[ji]
    var_lv3_south[:,ji] = var_lv3[t3,:,j,i]
    lon_lv3_south[ji] = lonr_lv3[j,i]


fig=plt.figure(figsize=(16,14))
gs = GridSpec(3,4)

#plot LV3 on left
# start with west boundary
vmin = Dvar['vmin']
vmax = Dvar['vmax']
ax0 = fig.add_subplot(gs[0,0])
p=ax0.pcolormesh(lat_lv3_west,sr_lv3,var_lv3_west,cmap='YlOrRd',vmin=vmin,vmax=vmax,shading='nearest')
cbaxes = inset_axes(ax0, width="4%", height="40%", loc=4,bbox_transform=ax0.transAxes,bbox_to_anchor=(-0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel(Dvar['axlabel'],rotation=90,labelpad=10,fontweight='bold')

ax0.set_ylabel('relative depth')
ax0.set_xlabel('latitude')
labeltext= 'LV3 \n'+date.strftime("%m/%d/%Y") + '\n'+ var_name + '\n' + 'west'
ax0.text(0.1,0.9,labeltext,transform=ax0.transAxes,fontweight='bold',va='top')
ax0.text(0.1,0.1,'shape = ({:},{:})'.format(var_lv3_west.shape[0],var_lv3_west.shape[1]),transform=ax0.transAxes,fontweight='bold',ha='left')

# next LV3 south boundary
ax1 = fig.add_subplot(gs[1,0])
ax1.pcolormesh(lon_lv3_south,sr_lv3,var_lv3_south,cmap='YlOrRd',vmin=vmin,vmax=vmax,shading='nearest')
ax1.set_ylabel('relative depth')
ax1.set_xlabel('longitude')
labeltext= 'LV3 \n'+date.strftime("%m/%d/%Y") + '\n'+ var_name + '\n' + 'south'
ax1.text(0.1,0.9,labeltext,transform=ax1.transAxes,fontweight='bold',va='top')
ax1.text(0.1,0.1,'shape = ({:},{:})'.format(var_lv3_south.shape[0],var_lv3_south.shape[1]),transform=ax1.transAxes,fontweight='bold',ha='left')

#finally add time series at bottom of western boundary
ax2 = fig.add_subplot(gs[2,:])
dt_list = []
nt3 = var_lv3.shape[0]
for tt in ott:
    dt_list.append(datetime(1999,1,1,0,0)+timedelta(seconds=tt))
j0,i0 = LV3_ji_west[0]
var_lv3_timeseries = var_lv3[:,-1,j0,i0]

ax2.plot(dt_list,var_lv3_timeseries,color='cornflowerblue',label='LV3')
ax2.plot(dt_list,var_west_LV4[t40:t40+nt3,-1,0],color='orange',label='LV4')
ax2.legend()
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax2.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax2.set_xlabel('Time')

# next LV4 west boundary
ax3 = fig.add_subplot(gs[0,1])
ax3.pcolormesh(latr_lv4[:,0],sr_lv4,var_west_LV4[t4,:,:],cmap='YlOrRd',vmin=vmin,vmax=vmax,shading='nearest')
ax3.set_ylabel('relative depth')
ax3.set_xlabel('latitude')
labeltext= 'LV4 \n'+date.strftime("%m/%d/%Y") + '\n'+ var_name + '\n' + 'west'
ax3.text(0.1,0.9,labeltext,transform=ax3.transAxes,fontweight='bold',va='top')
ax3.text(0.1,0.1,'shape = ({:},{:})'.format(var_west_LV4.shape[1],var_west_LV4.shape[2]),transform=ax3.transAxes,fontweight='bold',ha='left')

# next LV4 south boundary
ax4 = fig.add_subplot(gs[1,1])
ax4.pcolormesh(lonr_lv4[0,:],sr_lv4,var_south_LV4[t4,:,:],cmap='YlOrRd',vmin=vmin,vmax=vmax,shading='nearest')
ax4.set_ylabel('relative depth')
ax4.set_xlabel('longitude')
labeltext= 'LV4 \n'+date.strftime("%m/%d/%Y") + '\n'+ var_name + '\n' + 'south'
ax4.text(0.1,0.9,labeltext,transform=ax4.transAxes,fontweight='bold',va='top')
ax4.text(0.1,0.1,'shape = ({:},{:})'.format(var_south_LV4.shape[1],var_south_LV4.shape[2]),transform=ax4.transAxes,fontweight='bold',ha='left')

#compare horizontal resolution of LV3 and LV4
ax5 = fig.add_subplot(gs[0,2])
ax5.plot(lat_lv3_west,var_lv3_west[-1,:],color='cornflowerblue',label='LV3')
ax5.plot(lat_lv3_west,var_lv3_west[0,:],color='cornflowerblue',label='LV3')
ax5.plot(latr_lv4[:,0],var_west_LV4[t4,-1,:],color='orange',label='LV4')
ax5.plot(latr_lv4[:,0],var_west_LV4[t4,0,:],color='orange',label='LV4')
labeltext = var_name + ' @ surf & seafloor\nsnapshot\nwest'
ax5.text(0.1,0.9,labeltext,transform=ax5.transAxes,fontweight='bold',va='top')
ax5.set_xlabel('latitude')
ax5.set_ylabel(Dvar['axlabel'])
# ylim = ax5.get_ylim()
# ax5.text(0.1,0.9,'surface '+var_name+' snapshot along western boundary',transform=ax5.transAxes)

#compare horizontal resolution of LV3 and LV4
ax6 = fig.add_subplot(gs[1,2])
ax6.plot(lon_lv3_south,var_lv3_south[-1,:],color='cornflowerblue',label='LV3')
ax6.plot(lon_lv3_south,var_lv3_south[0,:],color='cornflowerblue',label='LV3')
ax6.plot(lonr_lv4[0,:],var_south_LV4[t4,-1,:],color='orange',label='LV4')
ax6.plot(lonr_lv4[0,:],var_south_LV4[t4,0,:],color='orange',label='LV4')
labeltext = var_name + ' @ surf & seafloor\nsnapshot\nsouth'
ax6.text(0.1,0.9,labeltext,transform=ax6.transAxes,fontweight='bold',va='top')
ax6.set_xlabel('longitude')
ax6.set_ylabel(Dvar['axlabel'])
ax6.set_ylim([Dvar['vmin'],Dvar['vmax']])
# ax6.text(0.1,0.9,'surface '+var_name+' snapshot along southern boundary',transform=ax6.transAxes)

# # compare values in LV3 and LV4
# ax7 = fig.add_subplot(gs[0,3])
# ax8 = fig.add_subplot(gs[1,3])
# for k in range(nz3):
#     ax7.scatter(lat_lv3_west,var_lv3_west[k,:],c='None',marker='o',edgecolors='cornflowerblue')
#     ax8.scatter(lon_lv3_south,var_lv3_south[k,:],c='None',marker='o',edgecolors='cornflowerblue')
# for k in range(nz4):
#     ax7.scatter(latr_lv4[:,0],var_west_LV4[t4,k,:],c='orange',marker='x')
#     ax8.scatter(lonr_lv4[0,:],var_south_LV4[t4,k,:],c='orange',marker='x')
#
# ax7.set_ylabel(var_name)
# ax7.set_xlabel('latitude')
# # ylim = ax7.get_ylim()
# ax8.set_ylabel(var_name)
# ax8.set_xlabel('longitude')
# ax8.set_ylim([16,17])

#add map to make sure you're extracting correctly
axmap = fig.add_subplot(gs[:2,-1])
#plot coastline
axmap.contour(lonr_lv3,latr_lv3,maskr_lv3,levels=[0.5],colors='k')
for ji in LV3_ji_west:
    j,i = ji
    axmap.scatter(lonr_lv3[j,i],latr_lv3[j,i],c='green')
for ji in LV3_ji_south:
    j,i = ji
    axmap.scatter(lonr_lv3[j,i],latr_lv3[j,i],c='magenta')
axmap.plot(lonr_lv3[j0,i0],latr_lv3[j0,i0],marker='*',markerfacecolor='orange',mec='k',markersize=15)
axmap.axis([-117.4,-117.05,32.3,32.8])
dar(axmap)
axmap.set_xlabel('longitude')
axmap.set_ylabel('latitude')

dlv3.close()
dlv4.close()
dgrd.close

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2019_BC_'+var_name+'_LV3_LV4_'+date.strftime("%Y.%m.%d")+'.png'
plt.savefig(out_fn)
