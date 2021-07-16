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
f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:14]=='ocean_his_NADB']

NOAA_tide_gauge_fname = '/data0/ebrasseale/WQ_data/2018validation/CO-OPS_9410170_met.csv'
NOAA_tide_gauge_df = pd.read_csv(NOAA_tide_gauge_fname,parse_dates={ 'time' : ['Date','Time (GMT)']})
NOAA_tide_gauge_df = NOAA_tide_gauge_df.set_index(NOAA_tide_gauge_df['time'])

NOAA_tide_gauge_lon = -117.17
NOAA_tide_gauge_lat = 32.71

NOAA_tide_gauge_ssh = NOAA_tide_gauge_df['Verified (m)'][:]
NOAA_tide_gauge_time = NOAA_tide_gauge_df['time']


ds = nc.Dataset(dir0+f_list[0])
lonr = ds['lon_rho'][:]
latr = ds['lat_rho'][:]
maskr = ds['mask_rho'][:]


latlondiff = np.sqrt((latr-NOAA_tide_gauge_lat)**2 + (lonr-NOAA_tide_gauge_lon+0.005)**2)
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

D = {}
var_list = ['CSIDE_time_list','CSIDE_ssh','NOAA_tide_gauge_ssh','NOAA_tide_gauge_time','NOAA_tide_gauge_lat','NOAA_tide_gauge_lon','iref','jref','lonr','latr','maskr']
for var_name in var_list:
    D[var_name] = locals()[var_name]


out_fn = '/data0/ebrasseale/WQ_data/CSIDE_2018_at_NOAA_tide_gauge.p'
pickle.dump(D,open(out_fn,'wb'))
