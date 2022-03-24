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
import pytz

parser = argparse.ArgumentParser()
parser.add_argument('-st', '--station', nargs='?', type=str, default='')
parser.add_argument('-yr', '--year', nargs='?', type=str, default='')
args = parser.parse_args()

# choose which year to extract
if len(args.station) == 0:
    print(30*'*' + ' x_CSIDE_at_dataset.py ' + 30*'*')
    print('\n%s\n' % '** Choose year (return for 2018) **')
    yr_list = ['2018', '2019']
    Nyr = len(yr_list)
    yr_dict = dict(zip(range(Nyr), yr_list))
    for nyr in range(Nyr):
        print(str(nyr) + ': ' + yr_list[nyr])
    my_nyr = input('-- Input number -- ')
    if len(my_nyr)==0:
        year = '2018'
    else:
        year = yr_dict[int(my_nyr)]
else:
    year = args.year

# choose which station to extract
if len(args.station) == 0:
    print(90*'*')
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


dir0 = '/data0/NADB'+year+'/'
f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[:14]=='ocean_his_NADB']
testing=False
if testing:
    f_list = f_list[:3]

if station=='NOAA':
    data_dict = {}
    data_dict['dataset_name'] = 'NOAA tide gauge 9410170 - San Diego, CA'
    # if year=='2018':
        # data_dict['fname'] = '/data0/ebrasseale/WQ_data/validation/CO-OPS_9410170_met_2018.csv'
    # elif year=='2019':
        # data_dict['fname'] = '/data0/ebrasseale/WQ_data/validation/CO-OPS_9410170_met_2019.csv'
    data_dict['fname'] = '/data0/ebrasseale/WQ_data/validation/CO-OPS_9410170_met_'+year+'.csv'
    data_dict['df'] = pd.read_csv(data_dict['fname'],parse_dates={ 'time' : ['Date','Time (GMT)']})
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['lon'] = -117.17
    data_dict['lat'] = 32.71
    data_dict['var_list'] = ['SSH (m)']
    var_list_df = {}
    var_list_df['SSH (m)'] = 'Verified (m)'
    var_list_roms = {}
    var_list_roms['SSH (m)'] = 'zeta'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = data_dict['df'][var_list_df[var_name]]

if station=='CDIP':
    data_dict = {}
    data_dict['dataset_name'] = '155 - Imperial Beach Nearshore Buoy (NDBC 46235)'
    data_dict['fname'] = '/data0/ebrasseale/WQ_data/validation/pm155p1p1_197501-202212.csv'
    dateparse = lambda x: datetime.strptime(x, '%Y %m %d %H %M')
    data_dict['df'] = pd.read_csv(data_dict['fname'],delim_whitespace=True,skiprows=[1],parse_dates={'time':['YEAR','MO','DY','HR','MN']},date_parser=dateparse)
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['lon'] = -117.16880
    data_dict['lat'] = 32.56957
    data_dict['var_list'] = ['Hs (m)','Tp (s)','Dp (deg)','SST (C)']
    var_list_df = {}
    var_list_df['Hs (m)'] = 'Hs'
    var_list_df['Tp (s)'] = 'Tp'
    var_list_df['Dp (deg)'] = 'Dp'
    var_list_df['SST (C)'] = 'Ta' # bc both delimiters and missing data are spaces, it's hard to parse across empty columns
    var_list_roms = {}
    var_list_roms['Hs (m)'] = 'Hwave'
    var_list_roms['Tp (s)'] = 'Pwave_top'
    var_list_roms['Dp (deg)'] = 'Dwave'
    var_list_roms['SST (C)'] = 'temp'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = data_dict['df'][var_list_df[var_name]]
        
if station=='SBOO':
    z_list = [1,10,18,26]
    ndepths = len(z_list)
    data_dict = {}
    data_dict['dataset_name'] = 'South Bay Ocean Outfall mooring'
    data_dict['fname_salt'] = '/data0/ebrasseale/WQ_data/validation/SBOO_sal_QC.csv'
    data_dict['fname_temp'] = '/data0/ebrasseale/WQ_data/validation/SBOO_temp_QC.csv'
    df_salt = pd.read_csv(data_dict['fname_salt'],parse_dates={  'time' : ['DateTime_PST']})
    df_temp = pd.read_csv(data_dict['fname_temp'],parse_dates={  'time' : ['DateTime_PST']})
    data_dict['df'] = pd.concat([df_salt, df_temp])
    data_dict['df'] = data_dict['df'].set_index(data_dict['df']['time'])
    data_dict['time'] = data_dict['df']['time']
    data_dict['time'] = data_dict['time'].tz_localize(pytz.timezone("America/Los_Angeles"),ambiguous=False).tz_convert(pytz.utc)
    data_dict['lon'] = -117.18612
    data_dict['lat'] = 32.53166
    data_dict['var_list'] = ['Temp (C)','Salt (psu)']
    var_list_df = {}
    var_list_df['Temp (C)'] = ['T_C_1m','T_C_10m','T_C_18m','T_C_26m']
    var_list_df['Salt (psu)'] = ['S_1m','S_10m','S_18m','S_26m']
    var_list_roms = {}
    var_list_roms['Temp (C)'] = 'temp'
    var_list_roms['Salt (psu)'] = 'salt'
    for var_name in data_dict['var_list']:
        data_dict[var_name] = {}
        for depth in range(ndepths):
            # first, nan out any bad data
            # find QC info for that variable/depth combo
            qname = 'Qual_'+var_list_df[var_name][depth]
            # set the variable/depth to nan everywhere the QC isn't 1
            data_dict['df'][var_list_df[var_name][depth]][data_dict['df'][qname]>1]=np.nan
            # read the resulting dataset into the dict entry for that depth
            data_dict[var_name][z_list[depth]] = data_dict['df'][var_list_df[var_name][depth]]
data_dict['station'] = station

ds = nc.Dataset(dir0+f_list[0])
var_list = ['lon_rho','lat_rho','mask_rho']
if station=='SBOO':
    var_list.extend(['h','s_rho','Cs_r','Vtransform','hc'])
for var in var_list:
    locals()[var] = ds[var][:]

#note add new modifier to keep form extracting on land for NOAA tide gauge - don't universally shift by 0.005!
latlondiff = np.sqrt((lat_rho-data_dict['lat'])**2 + (lon_rho-data_dict['lon'])**2)
#mask latlondiff before finding min
latlondiff[mask_rho==0] = np.nan
lld_nanmin = np.where(latlondiff==np.nanmin(latlondiff))
iref = lld_nanmin[1][0]
jref = lld_nanmin[0][0]


if station=='SBOO':
    #calculate vertical stretching array to help identify depths
    h0 = h[jref,iref]
    zr0 = (s_rho*hc + Cs_r*h0) / (hc + h0)
    nz = zr0.shape[0]
    zr0_rs = np.reshape(zr0,(1,nz))

NT = 0
CSIDE = {}
CSIDE['ot'] = np.array([])
for var_name in data_dict['var_list']:
    if station=='SBOO': # note: all variables for SBOO are defined at multiple depths
        CSIDE[var_name] = {}
        for z in z_list:
            CSIDE[var_name][z] = np.array([])
    else:
        CSIDE[var_name] = np.array([])

nf = len(f_list)
count=1
for fname in f_list:
    print(f'Working on file {count:d} of {nf:d}')
    ds = nc.Dataset(dir0+fname)

    ot = ds['ocean_time'][:]
    CSIDE['ot'] = np.append(CSIDE['ot'],ot)
    
    if station=='SBOO':
        #calculate z coordinates at mooring
        # adapated from Parker's zrfun.get_Z() using a priori info about the ROMS output I'm working with
        zeta0 = ds['zeta'][:] #always splice netcdf4 data after reading in
        zeta = zeta0[:,jref,iref]
        
        #nt is same for all but the last file. Calculate it once for the first N-1 files,
        # then calculate it again for the last file
        if fname==f_list[0] or f_list[-1]:
            nt = zeta.shape[0]
            zr0_tile = np.tile(zr0_rs,(nt,1))
        
        zeta_rs = np.reshape(zeta,(nt,1))
        zeta_tile = np.tile(zeta_rs,(1,nz))
        z_rho = zeta_tile + (zeta_tile + h0)*zr0_tile
        
        kref = np.zeros((nt,ndepths))
        for depth in range(ndepths):
            zref = z_list[depths] #note: these are positive, z_rho is negative
            kref[:,depth] = np.argmin(np.abs(z_rho+zref),axis=1)
    
    for var_name in data_dict['var_list']:
        var = ds[var_list_roms[var_name]][:]
        if len(var.shape)==3:
            #2d variable
            CSIDE[var_name] = np.append(CSIDE[var_name],var[:,jref,iref])
        elif len(var.shape)==4:
            #3d variable
            if station=='SBOO': # note all variables at SBOO are 3d
                for depth in ndepths:
                    CSIDE[var_name][z_list[depth]] = np.append(CSIDE[var_name][z_list[depth]],var[:,kref[:,depth],jref,iref])
            # in this case, use surface value only
            CSIDE[var_name] = np.append(CSIDE[var_name],var[:,-1,jref,iref])

    ds.close()
    count+=1


CSIDE['time'] = []
for t in CSIDE['ot']:
    date = datetime(1999,1,1)+timedelta(seconds=t)
    CSIDE['time'].append(date)

D = {}
# var_list = ['CSIDE_time_list','CSIDE_ssh','data_dict_ssh','data_dict_time','data_dict_lat','data_dict_lon','iref','jref','lonr','latr','maskr']
var_list = ['CSIDE','data_dict','lon_rho','lat_rho','mask_rho','iref','jref']
for var_name in var_list:
    D[var_name] = locals()[var_name]


out_fn = '/data0/ebrasseale/WQ_data/CSIDE_'+year+'_at_'+station+'.p'
pickle.dump(D,open(out_fn,'wb'))
