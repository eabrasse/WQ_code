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

rivfn = '/data0/ebrasseale/NADB2018/Input/river_tracer_4river_NADB2017_0.nc'

# rivfn = '/home/x1wu/SDTJRE_EPA/mfiles/Run2016_2017/LV4/River_input/NADB2017/river_tracer_4river_NADB2017_0.nc'

dsriv = nc.Dataset(rivfn)

rt = dsriv['river_time'][:]
rt_list = []
for t in rt:
    rt_list.append(datetime(1999,1,1)+timedelta(days=t))
    
# rQ = dsriv['river_transport'][:]
rT = dsriv['river_temp'][:]

nt,nriv = rT.shape

fig,ax=plt.subplots(4,1,sharex=True,figsize=(8,10))

ax[0].plot(rt_list,-rT[:,:5])
ax[0].text(0.9,0.9,'Tijuana River Estuary (1–5)',transform=ax[0].transAxes,ha='right')
ax[0].set_ylabel('temp (deg C)')
ax[0].get_xaxis().set_visible("false")

ax[1].plot(rt_list,-rT[:,5])
ax[1].set_ylabel('temp (deg C)')
ax[1].get_xaxis().set_visible("false")
ax[1].text(0.9,0.9,'Punta Bandera (6)',transform=ax[1].transAxes,ha='right')

ax[2].plot(rt_list,-rT[:,6:8])
ax[2].set_ylabel('temp (deg C)')
ax[2].get_xaxis().set_visible("false")
ax[2].text(0.9,0.9,'Sweetwater (7–8)',transform=ax[2].transAxes,ha='right')

ax[3].plot(rt_list,-rT[:,8])
ax[3].set_ylabel('temp (deg C)')
ax[3].get_xaxis().set_visible("false")
ax[3].text(0.9,0.9,'Otay (9)',transform=ax[3].transAxes,ha='right')

ax[-1].get_xaxis().set_visible("true")
ax[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
plt.setp( ax[-1].xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
ax[-1].set_xlabel('Time')

dsriv.close()

plt.tight_layout()

out_fn = '/data0/ebrasseale/WQ_plots/ROMS_rivertemp_2018.png'
# out_fn = '/data0/ebrasseale/WQ_plots/ROMS_rivertransport_2017.png'
plt.savefig(out_fn)
