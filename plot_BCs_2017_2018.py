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
# There is a working BC file for 2017,
# I want to compare it with the BC file generated for 2018
#
# Eventually I want to look at all different variables around
# the three open edges, but I hate looking at too many plots
# So I'm doing one variable x 3 directions x 2 time series at a time
    
plt.close('all')

BC_2018 = '/home/agilroy/Codes_XWu/Model_BCs/BC_LV4_20180101_20190101_Nz10.nc'
BC_2017 = '/home/x1wu/SDTJRE_EPA/mfiles/Run2016_2017/LV4/BC_LV4/BC_LV4_20170720_20180101_Nz10_dye.nc'

fig=plt.figure(figsize=(12,10))
gs = GridSpec(3,2)

var_name = 'ubar'

col = 0
for fn in BC_2017,BC_2018:
    ds = nc.Dataset(fn)
    tt = ds['temp_time'][:]

    dt_list = []
    for ott in tt:
        dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))
    
    row=0
    for edge in 'south','north','west':
        var = var_name+'_'+edge
        if row==0:
            ax0 = fig.add_subplot(gs[row,col])
            ax = ax0
            ax.title(['2017','2018'][col])
        else:
            ax = fig.add_subplot(gs[row,col],sharex=ax0)
            
        v = ds[var][:]
        y = np.arange(0,v.shape[1])
        p = ax.pcolormesh(dt_list,y,np.transpose(v),cmap='rainbow')
        ax.text(0.1,0.9,edge,color='black',fontweight='bold',transform=ax.transAxes)
        
        ax.set_ylabel('index along boundary')
        ax.get_xaxis().set_visible(False)
        
        if row==2:
            ax.get_xaxis().set_visible(True)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %Y"))
            plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
            ax.set_xlabel('Time',fontweight='bold')
        
        row+=1

    col+=1
    ds.close()

plt.tight_layout()
out_fn = '/data0/ebrasseale/WQ_plots/ROMS_2018_BC_'+var_name+'.png'
plt.savefig()
