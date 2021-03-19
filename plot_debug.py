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

rst_fn = home+'/NADB2018/ocean_rst_NADB2018_r1.nc'
his_fn = home+'/NADB2018/ocean_his_NADB2018_00007.nc'
riv_fn = home+'/NADB2018/Input/river_tracer_4river_NADB2018.nc'

dsr = nc.Dataset(rst_fn)
otdb = dsr['ocean_time'][:]
udb = dsr['u'][:]
ut = np.where(udb==udb.max())[0][0]
utwo = np.where(udb==udb.max())[1][0]
uz = np.where(udb==udb.max())[2][0]
uy = np.where(udb==udb.max())[3][0]
ux = np.where(udb==udb.max())[4][0]

ds = nc.Dataset(his_fn)
u = ds['u'][:]
lon_rho = ds['lon_rho'][:]
lat_rho = ds['lat_rho'][:]
mask_rho = ds['mask_rho'][:]

ot = ds['ocean_time'][:]
dt_list = []
for ott in ot:
    dt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))
    
dt_listdb = []
for ott in otdb:
    dt_listdb.append(datetime(1999,1,1,0,0) + timedelta(seconds=ott))

fig = plt.figure(figsize=(12,10))
gs = GridSpec(2,2)
axmap = fig.add_subplot(gs[:,0])
axmap.contour(lon_rho,lat_rho,mask_rho,colors='k',levels=[1],linewidths=0.5,alpha=1.0)
ucol='green'
axmap.plot(lon_rho[uy,ux],lat_rho[uy,ux],'*',mfc=ucol,mec='black',markersize=10)
axmap.text(lon_rho[uy,ux],lat_rho[uy,ux]+0.001,'max u',color=ucol,va='bottom',ha='center')

yl = axmap.get_ylim()
yav = (yl[0] + yl[1])/2
axmap.set_aspect(1/np.cos(np.pi*yav/180))
axmap.set_ylabel('Latitude')
axmap.set_xlabel('Longitude')

axu = fig.add_subplot(gs[0,1])
axu.plot(dt_list,u[:,uz,uy,ux])
axu.plot(dt_listdb,udb[:,utwo,uz,uy,ux],linestyle='dashed')
axu.set_title('u at location of blow up')
axu.set_xlabel('time')
axu.set_ylabel('velocity (m/s)')

axriv = fig.add_subplot(gs[1,1])
dsriv = nc.Dataset(riv_fn)
rt = dsriv['river_time'][:]
#convert from hours to seconds
rt = rt *3600*24
r0 = np.argmin(np.abs(rt-ot[0]))
r1 = np.argmin(np.abs(rt-ot[1]))
rt_list = []
for rtt in rt[r0:r1]:
    rt_list.append(datetime(1999,1,1,0,0) + timedelta(seconds=rtt))
rQ = dsriv['river_transport'][:]
axriv.plot(rt_list,rQ[r0:r1,:5])
axriv.text(0.9,0.9,'Tijuana River Estuary input',transform=axriv.transAxes,ha='right')
axriv.set_ylabel('transport m3s-1')
axriv.set_xlabel('Time')

plt.tight_layout()

outfn = home + 'WQ_plots/debug.png'
plt.savefig(outfn)
plt.close('all')

dsr.close()
dsriv.close()
ds.close()