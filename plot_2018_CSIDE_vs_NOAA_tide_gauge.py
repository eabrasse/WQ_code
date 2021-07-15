#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare CSIDE output with NOAA tide gauge
"""

# setup
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean as cmo
import netCDF4 as nc
import numpy as np
import os
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd
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

plt.close('all')
dir0 = '/data0/NADB2018/'
f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:14]=='ocean_his_NADB']

NOAA_tide_gauge_fname = '/data0/ebrasseale/WQ_data/2018validation/CO-OPS_9410170_met.csv'
NOAA_tide_gauge_df = pd.read_csv(NOAA_tide_gauge_fname,parse_dates={ 'time' : ['Date','Time (GMT)']})
NOAA_tide_gauge_df = NOAA_tide_gauge_df.set_index(NOAA_tide_gauge_df['time'])

NOAA_tide_gauge_lon = -117.17
NOAA_tide_gauge_lat = 32.71

NOAA_tide_gauge_ssh = NOAA_tide_gauge_df['Verified (m)'][:]


ds = nc.Dataset(dir0+f_list[0])
lonr = ds['lon_rho'][:]
latr = ds['lat_rho'][:]
maskr = ds['mask_rho'][:]


latlondiff = np.sqrt((latr-NOAA_tide_gauge_lat)**2 + (lonr-NOAA_tide_gauge_lon)**2)
iref = np.where(latlondiff==latlondiff.min())[1][0]
jref = np.where(latlondiff==latlondiff.min())[0][0]


NT = 0
CSIDE_time = np.array([])
CSIDE_ssh = np.array([])
for fname in f_list:
    ds = nc.Dataset(dir0+fname)
    nt = ds['ocean_time'].shape[0]
    NT +=nt
#     ds.close()
#
#
# old_nt = 0
# for fname in f_list:
#
#     ds = nc.Dataset(dir0+fname)
#
    ot = ds['ocean_time'][:]
    zeta = ds['zeta'][:]
    
    CSIDE_time = np.append(CSIDE_time,ot)
    CSIDE_ssh = np.append(CSIDE_ssh,zeta[:,jref,iref])
#
#     nt = ot.shape[0]
#
#     CSIDE_time[old_nt:old_nt+nt] = ot
#     CSIDE_ssh[old_nt:old_nt+nt] = zeta[:,jref,iref]
#
#     old_nt += nt
    ds.close()


CSIDE_time_list = []
for t in CSIDE_time:
    date = datetime(1999,1,1)+timedelta(seconds=t)
    CSIDE_time_list.append(date)


fig=plt.figure(figsize=(12,8))
gs = GridSpec(2,2)

# plot location of tide gauge
ax_map = fig.add_subplot(gs[:,1])
ax_map.contour(lonr,latr,maskr,levels=[0.5],colors='k')
ax_map.plot(lonr[jref,iref],latr[jref,iref],marker='*',markersize=15,mec='o',mfc='None')
ax_map.plot(lonr[jref,iref]-0.01,latr[jref,iref],marker='*',markersize=15,mec='k',mfc='yellow')
ax_map.set_xlabel('longitude')
ax_map.set_ylabel('latitude')
dar(ax_map)


# time series
ax_ts = fig.add_subplot(gs[0,0])
ax_ts.plot(CSIDE_time_list,CSIDE_ssh,label='CSIDE 2018')
ax_ts.plot(NOAA_tide_gauge_df['time'],NOAA_tide_gauge_ssh,label='NOAA tide gauge 9410170 - San Diego, CA')
ax_ts.legend()
ax_ts.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax_ts.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_ts.set_xlabel('Time')
ax_ts.set_ylabel('SSH (m)')
ax_ts.set_title('SSH at starred location')

# demeaned time series
ax_ts2 = fig.add_subplot(gs[1,0])
ax_ts2.plot(CSIDE_time_list,CSIDE_ssh-CSIDE_ssh.mean(),label='CSIDE 2018')
ax_ts2.plot(NOAA_tide_gauge_df['time'],NOAA_tide_gauge_ssh-NOAA_tide_gauge_ssh.mean(),label='NOAA tide gauge 9410170 - San Diego, CA')
# ax_ts2.legend()
ax_ts2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax_ts2.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_ts2.set_xlabel('Time')
ax_ts2.set_ylabel('SSH (m)')
ax_ts2.set_title('De-meaned SSH at starred location')

NOAA_tide_gauge_t0 = min(NOAA_tide_gauge_df['time'])
CSIDE_t0 = min(CSIDE_time_list)
later_t0 = max([NOAA_tide_gauge_t0,CSIDE_t0])

NOAA_tide_gauge_t1 = max(NOAA_tide_gauge_df['time'])
CSIDE_t1 = max(CSIDE_time_list)
earlier_t1 = min([NOAA_tide_gauge_t1,CSIDE_t1])

ssh1 = CSIDE_ssh[(CSIDE_time_list>later_t0)&(CSIDE_time_list<earlier_t1)]
ssh2 = NOAA_tide_gauge_ssh[(NOAA_tide_gauge_df['time']>later_t0)&(NOAA_tide_gauge_df['time']<earlier_t1)]
shared_time = CSIDE_time_list[(CSIDE_time_list>later_t0)&(CSIDE_time_list<earlier_t1)]
ssh_diff = ssh1-ssh2

ax_ts2_diff = ax_ts2.twinx()
ax_ts2_diff.plot(ssh_diff,color='k',linestyle='dotted',label='difference between demeaned SSHs')
ax_ts2_diff.legend()
ax_ts2_diff.set_ylabel(r'$\Delta$ SSH (m)')

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/CSIDE_2018_vs_NOAA_tide_gauge.png'
plt.savefig(out_fn)
