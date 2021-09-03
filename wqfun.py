#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experime.
"""

# setup
import os
import sys
alp = os.path.abspath('/Users/elizabethbrasseale/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cmocean as cmo
import pickle
import numpy as np

from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas
import matplotlib.transforms as mtrans
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import interp1d
import skill_metrics as sm


def willmott(m,o):
    """
    Calculates willmott skill score between two vectors, a model and a set of observations
    """
    # if len(m)!=len(o):
    #     error('Vectors must be the same length!');
    # end

    MSE = np.nanmean((m - o)**2)
    denom1 = abs(m - np.nanmean(o))
    denom2 = abs(o - np.nanmean(o))
    denom = np.nanmean((denom1 + denom2)**2)
    
    if denom==0:
        WS = 0
    else:
        WS = 1 - MSE/denom
    
    return WS

def earth_rad(lat_deg):
    """
    Calculate the Earth radius (m) at a latitude
    (from http://en.wikipedia.org/wiki/Earth_radius) for oblate spheroid

    INPUT: latitude in degrees

    OUTPUT: Earth radius (m) at that latitute
    """
    a = 6378.137 * 1000; # equatorial radius (m)
    b = 6356.7523 * 1000; # polar radius (m)
    cl = np.cos(np.pi*lat_deg/180)
    sl = np.sin(np.pi*lat_deg/180)
    RE = np.sqrt(((a*a*cl)**2 + (b*b*sl)**2) / ((a*cl)**2 + (b*sl)**2))
    return RE

def ll2xy(lon, lat, lon0, lat0):
    """
    This converts lon, lat into meters relative to lon0, lat0.
    It should work for lon, lat scalars or arrays.
    NOTE: lat and lon are in degrees!!
    """
    R = earth_rad(lat0)
    clat = np.cos(np.pi*lat0/180)
    x = R * clat * np.pi * (lon - lon0) / 180
    y = R * np.pi * (lat - lat0) / 180
    return x, y


def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.sin(np.pi*yav/180))
    
def get_shoreline_models(model_name_list):

    dir0 = '/Users/elizabethbrasseale/Projects/Water quality/WQ_data/'
    
    #load in data
    # #cside dye
    CSIDE_SAsm10 = {}
    cside_fn = dir0+'extractions2017/shoreline_dye_waves_05m_interp.p'
    Dcside = pickle.load(open(cside_fn,'rb'))
    CSIDE_SAsm10['dye'] = Dcside['dye_01'][:]
    CSIDE_SAsm10['y'] = Dcside['rshore'][:]

    #cside velocities
    cside_uv_fn = dir0 + 'extractions2017/shoreline_uv_05m_interp.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_SAsm10['v'] = Duv['v'][:]
    CSIDE_SAsm10['label']='CSIDE extracted data'# (rotated using\nsmoothed shoreangle (window=10))'

    CSIDE_PAsm100 = {}
    CSIDE_PAsm100['dye'] = Dcside['dye_01'][:]
    CSIDE_PAsm100['y'] = Dcside['rshore'][:]
    cside_uv_fn = dir0 + 'extractions2017/shoreline_uv_05m_interp_sm100_PA.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_PAsm100['v'] = Duv['v'][:]
    CSIDE_PAsm100['label']='CSIDE extracted data (rotated using\nsmoothed principle axis (window=100))'

    CSIDE_SAsm100 = {}
    CSIDE_SAsm100['dye'] = Dcside['dye_01'][:]
    CSIDE_SAsm100['y'] = Dcside['rshore'][:]
    cside_uv_fn = dir0 + 'extractions2017/shoreline_uv_05m_interp_sm100_SA.p'
    Duv = pickle.load(open(cside_uv_fn,'rb'))
    CSIDE_SAsm100['v'] = Duv['v'][:]
    CSIDE_SAsm100['label']='CSIDE extracted data'

    #adv-diff model wtih resolved alongshore-varying velocity w/ no smoothing
    AV_nosm = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_waves_postinterp.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_nosm['dye'] = da_model['c'][:]
    AV_nosm['y'] = da_model['y'][:]
    AV_nosm['v'] = da_model['v'][:]
    AV_nosm['label'] = 'Adv-Diff model with alongshore-varying input'#'\nusing smoothed shoreangle (window=10)'

    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_tuned['v'] = da_model['v'][:]
    AV_SAsm100_tuned['label'] = 'Adv-Diff model with alongshore-varying input'#'\nusing smoothed shoreangle (window=10)'

    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_3km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_3km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_3km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_3km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_3km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_3km_tuned['label'] = 'Adv-Diff model with alongshore-varying input, 3 km running average'#'\nusing smoothed shoreangle (window=10)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_4km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_4km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_4km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_4km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_4km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_4km_tuned['label'] = 'Adv-Diff model with alongshore-varying input, 4 km running average'#'\nusing smoothed shoreangle (window=10)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_5km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_5km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_5km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_5km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_5km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_5km_tuned['label'] = 'Adv-Diff model with alongshore-varying input, 5 km running average'#'\nusing smoothed shoreangle (window=10)'
    
    #adv-diff model wtih resolved alongshore-varying velocity, rotated using SA sm100 and tuned
    AV_SAsm100_6km_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA_6km_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100_6km_tuned['dye'] = da_model['c'][:]
    AV_SAsm100_6km_tuned['y'] = da_model['y'][1:-1]
    AV_SAsm100_6km_tuned['v'] = da_model['v'][:]
    AV_SAsm100_6km_tuned['label'] = 'Adv-Diff model with alongshore-varying input, 6 km running average'#'\nusing smoothed shoreangle (window=10)'

    #adv-diff model with CSIDE-extracted velocity input
    AV_recycled_tuned = {}
    model_fn=dir0+'adv_diff_model/CSIDE_recycled_input_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_recycled_tuned['y'] = da_model['y'][:]
    AV_recycled_tuned['dye']  = da_model['c'][:]
    AV_recycled_tuned['v']  = da_model['v'][:]
    AV_recycled_tuned['label'] = 'Adv-Diff model with velocity directly from CSIDE'

    #adv-diff model with 3km binned velocity
    AV_SAsm10_3kmbin = {}
    model_fn=dir0+'adv_diff_model/CSIDE_3kmbin_waves_postinterp.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_3kmbin['y'] = da_model['y'][:]
    AV_SAsm10_3kmbin['dye']  = da_model['c'][:]
    AV_SAsm10_3kmbin['v']  = da_model['v'][:]
    AV_SAsm10_3kmbin['label'] = 'Adv-Diff model with 3km-binned velocity input'

    #adv-diff model with 1km running average velocity
    AV_SAsm10_1kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_1km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_1kmavg['y'] = da_model['y'][:]
    AV_SAsm10_1kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_1kmavg['v']  = da_model['v'][:]
    AV_SAsm10_1kmavg['label'] = 'Adv-Diff model with 1km-running avg velocity input'

    #adv-diff model with 3km running average velocity
    AV_SAsm10_3kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_3km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_3kmavg['y'] = da_model['y'][:]
    AV_SAsm10_3kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_3kmavg['v']  = da_model['v'][:]
    AV_SAsm10_3kmavg['label'] = 'Adv-Diff model with 3km-running avg velocity input'

    #adv-diff model with 5km running average velocity
    AV_SAsm10_5kmavg = {}
    model_fn=dir0+'adv_diff_model/CSIDE_runavg_5km.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm10_5kmavg['y'] = da_model['y'][:]
    AV_SAsm10_5kmavg['dye']  = da_model['c'][:]
    AV_SAsm10_5kmavg['v']  = da_model['v'][:]
    AV_SAsm10_5kmavg['label'] = 'Adv-Diff model with 5km-running avg velocity input'


    #adv-diff model with uniform velocity
    U_tuned = {}
    model_fn2=dir0+'adv_diff_model/CSIDE_tuning_kd_-2.67E-05.p'
    da_model2 = pickle.load(open(model_fn2,'rb'))
    U_tuned['y']= da_model2['y'][:]
    U_tuned['dye'] = da_model2['c'][:]
    v00 = da_model2['v'][:]
    nt,nj = U_tuned['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_tuned['v'] = v0
    U_tuned['label'] = 'Adv-Diff model with alongshore uniform velocity input'

    #adv-diff model wtih resolved alongshore-varying velocity
    #but leaving waves uninterpolated, and interpolating v after calculating Sxy
    AV_PAsm100 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_PA.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_PAsm100['dye'] = da_model['c'][:]
    AV_PAsm100['y'] = da_model['y'][1:-1]
    AV_PAsm100['v'] = da_model['v'][:]
    AV_PAsm100['label'] = 'Adv-Diff model with alongshore-varying input\nusing smoothed PA (window=100)'

    #adv-diff model wtih resolved alongshore-varying velocity
    #but leaving waves uninterpolated, and interpolating v after calculating Sxy
    AV_SAsm100 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_alongshore_varying_sm100_SA.p'
    da_model = pickle.load(open(model_fn,'rb'))
    AV_SAsm100['dye'] = da_model['c'][:]
    AV_SAsm100['y'] = da_model['y'][1:-1]
    AV_SAsm100['v'] = da_model['v'][:]
    AV_SAsm100['label'] = 'Adv-Diff model with alongshore-varying input\nusing smoothed shoreangles (window=100)'


    #alongshore diffusion model
    U_Kyy01 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy01_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy01['dye'] = da_model['c'][:,1:-1]
    U_Kyy01['y'] = da_model['y'][1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_Kyy01['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy01['v'] = v0
    U_Kyy01['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 1 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy03 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy03_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy03['dye'] = da_model['c'][:,1:-1]
    U_Kyy03['y'] = da_model['y'][1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_Kyy03['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy03['v'] = v0
    U_Kyy03['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 3 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy05 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy05_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy05['dye'] = da_model['c'][:,1:-1]
    U_Kyy05['y'] = da_model['y'][1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_Kyy05['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy05['v'] = v0
    U_Kyy05['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 5 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy07 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy07_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy07['dye'] = da_model['c'][:,1:-1]
    U_Kyy07['y'] = da_model['y'][1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_Kyy07['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy07['v'] = v0
    U_Kyy07['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 7 m^{2}s^{-1}$'
    
    #alongshore diffusion model
    U_Kyy10 = {}
    model_fn=dir0+'adv_diff_model/CSIDE_PB_Kyy10_tuned.p'
    da_model = pickle.load(open(model_fn,'rb'))
    U_Kyy10['dye'] = da_model['c'][:,1:-1]
    U_Kyy10['y'] = da_model['y'][1:-1]
    v00 = da_model2['v'][:]
    nt,nj = U_Kyy10['dye'].shape
    v00 = np.reshape(v00,(nt,1))
    v0 = np.tile(v00,(1,nj))
    U_Kyy10['v'] = v0
    U_Kyy10['label'] = r'Adv-Diff model alongshore diffusion $Kyy = 10 m^{2}s^{-1}$'

    model_dict_list = {}

    for var in model_name_list:
        model_dict_list[var] = locals()[var]
    
    return(model_dict_list)