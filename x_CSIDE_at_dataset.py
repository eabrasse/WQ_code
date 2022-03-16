#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare CSIDE output with NOAA tide gauge
"""

# setup
import netCDF4 as nc
import numpy as np
import os
import argparse
from datetime import datetime, timedelta
import pandas as pd
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-st', '--station', nargs='?', type=str, default='')
args = parser.parse_args()

# choose which station to extract
if len(args.station) == 0:
    print(30*'*' + ' x_CSIDE_at_dataset.py ' + 30*'*')
    print('\n%s\n' % '** Choose station (return for NOAA) **')
    st_list = ['NOAA', 'CDIP', 'SBOO']
    Nst = len(st_list)
    st_dict = dict(zip(range(Nst), st_list))
    for nst in range(Nst):
        print(str(nst) + ': ' + st_list[nst])
    my_nst = input('-- Input number -- ')
    if len(my_nst)==0:
        station = 'NOAA'
    else:
        station = st_dict[int(my_nst)]
else:
    station = args.station

dir0 = '/data0/NADB2018/'
f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:14]=='ocean_his_NADB']


if station=='NOAA':
    data_dict = {}
    data_dict['dataset_name'] = 'NOAA tide gauge 9410170 - San Diego, CA'
    data_dict['fname'] = '/data0/ebrasseale/WQ_data/2018validation/CO-OPS_9410170_met.csv'
    data_dict['df'] = pd.read_csv(data_dict['fname'],parse_dates={ 'time' : ['Date','Time (GMT)']})
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['lon'] = -117.17
    data_dict['lat'] = 32.71
    data_dict['var_list'] = 'ssh'
    var_list_df = {}
    var_list_df['SSH'] = 'Verified (m)'
    var_list_roms = {}
    var_list_roms['SSH'] = 'zeta'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = data_dict['df'][var_list_df[var_name]]

if station=='CDIP':
    data_dict = {}
    data_dict['dataset_name'] = 'NOAA tide gauge 9410170 - San Diego, CA'
    data_dict['fname'] = '/data0/ebrasseale/WQ_data/2018validation/pm191p1p1_201801-201812_new.csv'
    data_dict['df'] = pd.read_csv(data_dict['fname'],parse_dates={ 'time' : ['Date','Time (GMT)']})
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['lon'] = -117.17
    data_dict['lat'] = 32.71
    data_dict['var_list'] = 'Hs','Tp','Dp','temp'
    var_list_df = {}
    var_list_df['Hs'] = 'Hs (m)'
    var_list_df['Tp'] = 'Tp (sec)'
    var_list_df['Dp'] = 'Dp (deg)'
    var_list_df['temp'] = 'Temp (Sfc (C))'
    var_list_roms = {}
    var_list_roms['Hs'] = 'Hwave'
    var_list_roms['Tp'] = 'Pwave_top'
    var_list_roms['Dp'] = 'Dwave'
    var_list_roms['temp'] = 'temp'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = data_dict['df'][var_list_df[var_name]]
        
if station=='SBOO':
    data_dict = {}
    data_dict['dataset_name'] = 'South Bay Ocean Outfall mooring'
    data_dict['fname_salt'] = '/data0/ebrasseale/WQ_data/2018validation/SBOO_sal_QC.csv'
    data_dict['fname_temp'] = '/data0/ebrasseale/WQ_data/2018validation/SBOO_temp_QC.csv'
    df_salt = pd.read_csv(data_dict['fname_salt'],parse_dates={  time : ['DateTime_PST']})
    df_temp = pd.read_csv(data_dict['fname_temp'],parse_dates={  time : ['DateTime_PST']})
    data_dict['df'] = pd.concat([df_salt, df_temp])
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['lon'] = -117.18612
    data_dict['lat'] = 32.53166
    data_dict['var_list'] = 'SST','SSS'
    var_list_df = {}
    var_list_df['SST'] = 'T_C_1m'
    var_list_df['SSS'] = 'S_1m'
    var_list_roms = {}
    var_list_roms['SST'] = 'temp'
    var_list_roms['SSS'] = 'salt'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = data_dict['df'][var_list_df[var_name]]

ds = nc.Dataset(dir0+f_list[0])
lonr = ds['lon_rho'][:]
latr = ds['lat_rho'][:]
maskr = ds['mask_rho'][:]


latlondiff = np.sqrt((latr-data_dict['lat'])**2 + (lonr-data_dict['lon']+0.005)**2)
iref = np.where(latlondiff==latlondiff.min())[1][0]
jref = np.where(latlondiff==latlondiff.min())[0][0]


NT = 0
CSIDE = {}
CSIDE['time'] = np.array([])
for var_name in data_dict['var_list']:
    CSIDE[var_name] = np.array([])
    
for fname in f_list:
    ds = nc.Dataset(dir0+fname)

    ot = ds['ocean_time'][:]
    CSIDE['time'] = np.append(CSIDE['time'],ot)
    
    for var_name in data_dict['var_list']:
        var = ds[var_list_roms[var_name]][:]
        CSIDE[var_name] = np.append(CSIDE[var_name],var[:,jref,iref])

    ds.close()


CSIDE['time_list'] = []
for t in CSIDE['time']:
    date = datetime(1999,1,1)+timedelta(seconds=t)
    CSIDE['time_list'].append(date)

D = {}
# var_list = ['CSIDE_time_list','CSIDE_ssh','data_dict_ssh','data_dict_time','data_dict_lat','data_dict_lon','iref','jref','lonr','latr','maskr']
var_list = ['CSIDE','data_dict','lonr','latr','maskr','iref','jref']
for var_name in var_list:
    D[var_name] = locals()[var_name]


out_fn = '/data0/ebrasseale/WQ_data/CSIDE_2018_at_'+station+'.p'
pickle.dump(D,open(out_fn,'wb'))
