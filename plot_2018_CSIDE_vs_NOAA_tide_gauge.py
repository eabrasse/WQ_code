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
import pickle

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))

plt.close('all')
dir0 = '/data0/NADB2018/'

data_fn = '/data0/ebrasseale/WQ_data/CSIDE_2018_at_NOAA_tide_gauge.p'
D = pickle.load(open(data_fn,'rb'))

for var_name in D.keys():
    locals()[var_name] = D[var_name]

fig=plt.figure(figsize=(12,8))
gs = GridSpec(2,2)

# plot location of tide gauge
ax_map = fig.add_subplot(gs[:,1])
ax_map.contour(lonr,latr,maskr,levels=[0.5],colors='k')
ax_map.plot(lonr[jref,iref],latr[jref,iref],marker='o',markersize=15,mec='k',mfc='None')
ax_map.plot(lonr[jref,iref]-0.01,latr[jref,iref],marker='*',markersize=15,mec='k',mfc='yellow')
ax_map.set_xlabel('longitude')
ax_map.set_ylabel('latitude')
dar(ax_map)


# time series
ax_ts = fig.add_subplot(gs[0,0])
ax_ts.plot(CSIDE_time_list,CSIDE_ssh,label='CSIDE 2018')
ax_ts.plot(NOAA_tide_gauge_time,NOAA_tide_gauge_ssh,label='NOAA tide gauge 9410170 - San Diego, CA')
ax_ts.legend()
ax_ts.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax_ts.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_ts.set_xlabel('Time')
ax_ts.set_ylabel('SSH (m)')
ax_ts.set_title('SSH at starred location')

# demeaned time series
ax_ts2 = fig.add_subplot(gs[1,0])
ax_ts2.plot(CSIDE_time_list,CSIDE_ssh-CSIDE_ssh.mean(),label='CSIDE 2018')
ax_ts2.plot(NOAA_tide_gauge_time,NOAA_tide_gauge_ssh-NOAA_tide_gauge_ssh.mean(),label='NOAA tide gauge 9410170 - San Diego, CA')
# ax_ts2.legend()
ax_ts2.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax_ts2.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax_ts2.set_xlabel('Time')
ax_ts2.set_ylabel('SSH (m)')
ax_ts2.set_title('De-meaned SSH at starred location')

NOAA_tide_gauge_t0 = min(NOAA_tide_gauge_time)
CSIDE_t0 = min(CSIDE_time_list)
t0 = max([NOAA_tide_gauge_t0,CSIDE_t0])

NOAA_tide_gauge_t1 = max(NOAA_tide_gauge_time)
CSIDE_t1 = max(CSIDE_time_list)
t1 = min([NOAA_tide_gauge_t1,CSIDE_t1])

# if type(t0) is pd.Timestamp:
#     t0 = t0.to_pydatetime()
# if type(t1) is pd.Timestamp:
#     t1 = t1.to_pydatetime()
shared_time = [x for x in CSIDE_time_list if x>t0 and x<t1]

ssh1 = np.array([CSIDE_ssh[i] for i in range(len(CSIDE_time_list)) if CSIDE_time_list[i]>t0 and CSIDE_time_list[i]<t1])
ssh2 = np.array([NOAA_tide_gauge_ssh[i] for i in range(len(NOAA_tide_gauge_time)) if NOAA_tide_gauge_time[i]>t0 and NOAA_tide_gauge_time[i]<t1])

ssh_diff = ssh1-ssh2
#
# ax_ts2_diff = ax_ts2.twinx()
ax_ts2.plot(shared_time,ssh_diff,color='k',linestyle='dotted',label='difference between demeaned SSHs')
ax_ts2.legend()
# ax_ts2_diff.set_ylabel(r'$\Delta$ SSH (m)')

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/CSIDE_2018_vs_NOAA_tide_gauge.png'
plt.savefig(out_fn)
