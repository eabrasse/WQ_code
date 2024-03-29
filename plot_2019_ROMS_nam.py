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
import wqfun
    
# Build a code to plot BC's generated by code
#
# I want to compare it with the BC file generated for 2018
# with the LV3 files for 2018
    
plt.close('all')
dir0 = '/data0/NADB2019/Input/'
LV4_grid = dir0+'GRID_SDTJRE_LV4_ROTATE_rx020_hplus020_DK_4river_otaymk.nc'
dir0_nam = '/home/x1wu/SDTJRE_2018/mfiles/NAM_data/LV3_2020/'

ts0 = {}
ts0['LV4_nam_2019'] = dir0+'roms_nam_LV4_20181216_20190101.nc'
ts0['LV3_nam_2019'] = dir0_nam+ 'roms_nam_LV3_20181216_20190115.nc'
# now match time indexes between the two grids
ts0['date0'] = datetime(2018,12,16)
ts0['date1'] = datetime(2019,1,1)


ts1 = {}
ts1['LV4_nam_2019'] = dir0+'roms_nam_LV4_20181231_20190120.nc'
ts1['LV3_nam_2019'] = dir0_nam+ 'roms_nam_LV3_20181216_20190115.nc'
ts1['date0'] = datetime(2018,12,31)
ts1['date1'] = datetime(2019,1,15)

ts2 = {}
ts2['LV4_nam_2019'] = dir0+'roms_nam_LV4_20190430_20190520.nc'
ts2['LV3_nam_2019'] = dir0_nam+ 'roms_nam_LV3_20190415_20190515.nc'
ts2['date0'] = datetime(2019,4,30)
ts2['date1'] = datetime(2019,5,15)

ts3 = {}
ts3['LV4_nam_2019'] = dir0+'roms_nam_LV4_20191206_20191226.nc'
ts3['LV3_nam_2019'] = dir0_nam+ 'roms_nam_LV3_20191211_20200101.nc'
# now match time indexes between the two grids
ts3['date0'] = datetime(2019,12,11)
ts3['date1'] = datetime(2019,12,26)

ts_list = ts0,ts1,ts2,ts3
count = 0
print('Choose time series to plot based on start time (return=0)')
for ts in ts_list:
    print('{}. '.format(count)+ts['date0'].strftime("%Y.%m.%d"))
    count+=1
user_input = input()
if len(user_input)==0:
    ts_ind = 0
else:
    ts_ind = int(user_input)
my_ts = ts_list[ts_ind]

dgrd = nc.Dataset(LV4_grid)
lonr_LV4 = dgrd['lon_rho'][:]
latr_LV4 = dgrd['lat_rho'][:]
maskr_LV4 = dgrd['mask_rho'][:]

dlv4 = nc.Dataset(my_ts['LV4_nam_2019'])
uwind_LV4 = dlv4['Uwind'][:]
vwind_LV4 = dlv4['Vwind'][:]
pair_LV4 = dlv4['Pair'][:]

dlv3 = nc.Dataset(my_ts['LV3_nam_2019'])
lonr_LV3 = dlv3['lon'][:]
latr_LV3 = dlv3['lat'][:]
uwind_LV3 = dlv3['Uwind'][:]
vwind_LV3 = dlv3['Vwind'][:]
pair_LV3 = dlv3['Pair'][:]

# # choose a spatial reference index
lonref=-117.2
latref = 32.55
latlondiff_LV3 = np.sqrt((latr_LV3-latref)**2 + (lonr_LV3-lonref)**2)
i3 = np.where(latlondiff_LV3==latlondiff_LV3.min())[1][0]
j3 = np.where(latlondiff_LV3==latlondiff_LV3.min())[0][0]

latlondiff_LV4 = np.sqrt((latr_LV4-latref)**2 + (lonr_LV4-lonref)**2)
i4 = np.where(latlondiff_LV4==latlondiff_LV4.min())[1][0]
j4 = np.where(latlondiff_LV4==latlondiff_LV4.min())[0][0]

wt3 = dlv3['wind_time'][:]
wt3_list = []
for t in wt3:
    date = datetime(1999,1,1)+timedelta(days=t)
    wt3_list.append(date)
t30 = wt3_list.index(my_ts['date0'])
t31 = wt3_list.index(my_ts['date1'])

wt4 = dlv4['wind_time'][:]
wt4_list = []
for t in wt4:
    date = datetime(1999,1,1)+timedelta(days=t)
    wt4_list.append(date)
t40 = wt4_list.index(my_ts['date0'])
t41 = wt4_list.index(my_ts['date1'])

fig=plt.figure(figsize=(12,16))
gs = GridSpec(3,3)

# #plot LV3 on left
# # plot u wind
vmin = -5
vmax = 5
ax0 = fig.add_subplot(gs[0,0])
p=ax0.pcolormesh(lonr_LV3,latr_LV3,uwind_LV3[t30,:],cmap='BrBG',vmin=vmin,vmax=vmax,shading='auto')
cbaxes = inset_axes(ax0, width="4%", height="40%", loc=4,bbox_transform=ax0.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('wind (ms-1)',rotation=90,labelpad=10,fontweight='bold')

ax0.set_ylabel('latitude')
ax0.set_xlabel('longitude')
labeltext= 'LV3 \n'+wt3_list[t30].strftime("%m/%d/%Y") + '\nUwind'
ax0.text(0.1,0.1,labeltext,transform=ax0.transAxes,fontweight='bold',va='top')

# plot LV3 v wind
ax1 = fig.add_subplot(gs[1,0])
p=ax1.pcolormesh(lonr_LV3,latr_LV3,vwind_LV3[t30,:],cmap='BrBG',vmin=vmin,vmax=vmax,shading='auto')
ax1.set_ylabel('latitude')
ax1.set_xlabel('longitude')
labeltext= 'LV3 \n'+wt3_list[t30].strftime("%m/%d/%Y") + '\nVwind'
ax1.text(0.1,0.1,labeltext,transform=ax1.transAxes,fontweight='bold',va='top')
cbaxes = inset_axes(ax1, width="4%", height="40%", loc=4,bbox_transform=ax1.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('wind (ms-1)',rotation=90,labelpad=10,fontweight='bold')

# plot LV3 air pressure
ax2 = fig.add_subplot(gs[2,0])
p2=ax2.pcolormesh(lonr_LV3,latr_LV3,pair_LV3[t30,:],cmap='RdPu',vmin=1013,vmax=1021,shading='auto')
ax2.set_ylabel('latitude')
ax2.set_xlabel('longitude')
labeltext= 'LV3 \n'+wt3_list[t30].strftime("%m/%d/%Y") + '\nPair'
ax2.text(0.1,0.1,labeltext,transform=ax2.transAxes,fontweight='bold',va='top')
cbaxes = inset_axes(ax2, width="4%", height="40%", loc=4,bbox_transform=ax2.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p2, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('pressure (Pa)',rotation=90,labelpad=10,fontweight='bold')

# LV4 uwind
ax3 = fig.add_subplot(gs[0,1])
p=ax3.pcolormesh(lonr_LV4,latr_LV4,uwind_LV4[t40,:],cmap='BrBG',vmin=vmin,vmax=vmax,shading='auto')
ax3.set_ylabel('latitude')
ax3.set_xlabel('longitude')
labeltext= 'LV4 \n'+wt4_list[t40].strftime("%m/%d/%Y") + '\nUwind'
ax3.text(0.1,0.1,labeltext,transform=ax3.transAxes,fontweight='bold',va='top')
cbaxes = inset_axes(ax3, width="4%", height="40%", loc=4,bbox_transform=ax3.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('wind (ms-1)',rotation=90,labelpad=10,fontweight='bold')

# LV4 vwind
ax4 = fig.add_subplot(gs[1,1])
p=ax4.pcolormesh(lonr_LV4,latr_LV4,vwind_LV4[t40,:,:],cmap='BrBG',vmin=vmin,vmax=vmax,shading='auto')
ax4.set_ylabel('latitude')
ax4.set_xlabel('longitude')
labeltext= 'LV4 \n'+wt4_list[t40].strftime("%m/%d/%Y") + '\nVwind'
ax4.text(0.1,0.1,labeltext,transform=ax4.transAxes,fontweight='bold',va='top')
cbaxes = inset_axes(ax4, width="4%", height="40%", loc=4,bbox_transform=ax4.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('wind (ms-1)',rotation=90,labelpad=10,fontweight='bold')

# LV4 air pressure
ax5 = fig.add_subplot(gs[2,1])
p=ax5.pcolormesh(lonr_LV4,latr_LV4,pair_LV4[t40,:],cmap='RdPu',vmin=1013,vmax=1021,shading='auto')
ax5.set_ylabel('latitude')
ax5.set_xlabel('longitude')
labeltext= 'LV4 \n'+wt4_list[t40].strftime("%m/%d/%Y") + '\nPair'
ax5.text(0.1,0.1,labeltext,transform=ax5.transAxes,fontweight='bold',va='top')
cbaxes = inset_axes(ax5, width="4%", height="40%", loc=4,bbox_transform=ax5.transAxes,bbox_to_anchor=(0.15,0.0,1,1))
cb = fig.colorbar(p, cax=cbaxes, orientation='vertical')
cb.ax.set_ylabel('pressure (Pa)',rotation=90,labelpad=10,fontweight='bold')

for ax in ax0,ax1,ax2,ax3,ax4,ax5:
    # for ax in ax3,ax4,ax5:
    #zoom in on LV4 axis limits
    ax.axis([lonr_LV4.min(),lonr_LV4.max(),latr_LV4.min(),latr_LV4.max()])
    
    #add the spatial reference point for the line plots
    ax.plot(lonr_LV4[j4,i4],latr_LV4[j4,i4],marker='*',markersize=10,markerfacecolor='yellow',markeredgecolor='k')
    
    #contour the shoreline and adjust aspect ratio
    ax.contour(lonr_LV4,latr_LV4,maskr_LV4,levels=[0.5],colors='k')
    wqfun.dar(ax)

# uwind time compare
ax6 = fig.add_subplot(gs[0,2])
ax6.plot(wt3_list[t30:t31],uwind_LV3[t30:t31,j3,i3],color='cornflowerblue',label='LV3')
ax6.plot(wt4_list[t40:t41],uwind_LV4[t40:t41,j4,i4],color='orange',label='LV4')
ax6.legend()
ax6.text(0.1,0.1,'Uwind time series @ star',transform=ax6.transAxes)
ax6.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax6.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax6.set_xlabel('Time')
ax6.set_ylabel('wind velocity (m/s)')

# vwind time compare
ax7 = fig.add_subplot(gs[1,2])
ax7.plot(wt3_list[t30:t31],vwind_LV3[t30:t31,j3,i3],color='cornflowerblue',label='LV3')
ax7.plot(wt4_list[t40:t41],vwind_LV4[t40:t41,j4,i4],color='orange',label='LV4')
ax7.legend()
ax7.text(0.1,0.1,'Vwind time series @ star',transform=ax7.transAxes)
ax7.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax7.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax7.set_xlabel('Time')
ax7.set_ylabel('wind velocity (m/s)')


# pair time compare
ax8 = fig.add_subplot(gs[2,2])
ax8.plot(wt3_list[t30:t31],pair_LV3[t30:t31,j3,i3],color='cornflowerblue',label='LV3')
ax8.plot(wt4_list[t40:t41],pair_LV4[t40:t41,j4,i4],color='orange',label='LV4')
ax8.legend()
ax8.text(0.1,0.1,'Pair time series @ star',transform=ax8.transAxes)
ax8.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax8.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax8.set_xlabel('Time')
ax8.set_ylabel('air pressure (Pa)')

dlv3.close()
dlv4.close()
dgrd.close

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2019_NAM_LV3_LV4_'+my_ts['date0'].strftime("%Y.%m.%d")+'.png'
# out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2018_NAM_LV4_'+date0.strftime("%Y.%m.%d")+'.png'
plt.savefig(out_fn)
