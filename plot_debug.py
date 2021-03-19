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
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta

home = '/data0/ebrasseale/'

fn = home+'/NADB2018/ocean_his_NADB2018_00007.nc'

ds = nc.Dataset(fn)
u = ds['u'][:]
ut = np.where(u==u.max)[0][0]
uz = np.where(u==u.max)[1][0]
uy = np.where(u==u.max)[2][0]
ux = np.where(u==u.max)[3][0]

zeta = ds['zeta'][:]
zetat = np.where(u==u.max)[0][0]
zetaz = np.where(u==u.max)[1][0]
zetay = np.where(u==u.max)[2][0]
zetax = np.where(u==u.max)[3][0]

lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
mask_rho = ds['mask_rho'][:]

ot = ds['ocean_time'][:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

fig = plt.figure(figsize=(12,10))
gs = GridSpec(2,2)
axmap = fig.add_subplots(gs[0,:])
axmap.contour(lon_rho,lat_rho,mask_rho,colors='k',levels=[1],linewidths=0.5,alpha=1.0)
ucol='green'
axmap.plot(lon_rho[uy,ux],lat_rho[uy,ux],'*',mfc=ucol,mec='black',markersize=10)
axmap.text(lon_rho[uy,ux],lat_rho[uy,ux]+0.001,'max u',color=ucol,va='bottom',ha='center')
zetacol='magenta'
axmap.plot(lon_rho[zetay,zetax],lat_rho[zetay,zetax],'d',mfc=zetacol,mec='black',markersize=10)
axmap.text(lon_rho[zetay,zetax],lat_rho[zetay,zetax]-0.001,'max SSH',color=zetacol,va='top',ha='center')

yl = axmap.get_ylim()
yav = (yl[0] + yl[1])/2
axmap.set_aspect(1/np.cos(np.pi*yav/180))
axmap.set_ylabel('Latitude')
axmap.set_xlabel('Longitude')

axu = fig.add_subplots(gs[1,0])
axu.plot(dt_list,u[:,uz,uy,ux])
axu.set_title('u at location of blow up')
axu.set_xlabel('time')
axu.set_ylabel('velocity (m/s)')

axzeta = fig.add_subplots(gs[1,1])
axzeta.plot(dt_list,zeta[:,zetay,zetax])
axzeta.set_title('zeta at location of blow up')
axzeta.set_xlabel('time')
axzeta.set_ylaebl('SSH (m)')

plt.tight_layout()

outfn = home + 'WQ_plots/debug.png'
plt.savefig(outfn)
plt.close('all')
