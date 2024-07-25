import io
import re
import sys
import time
import requests
import os
import pandas as pd
import numpy as np
from io import StringIO
import matplotlib.pyplot as plt

from lightcurvy import util
from mag_to_flux import mag_to_flux

###

def panstarrs_query(obj = None, ra = None, dec = None, search_radius = 0.0015, page_size = 100,
                    save = True, overwrite = False, obj_save = None, save_dir = None, **kwargs):
    
    print('Grabbing Pan-STARRS light curve')
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'panstarrs', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'panstarrs'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'Pan-STARRS light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        url = rf'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv?ra={ra}&dec={dec}&radius={search_radius}&pagesize={page_size}&format=csv'
        
        print(f'\tObject: {obj}\n\tSearch Radius: {search_radius}\n\tURL: {url}')
        
        response = requests.get(url)
        
        try:
            lc = pd.read_csv(StringIO(response.text))
        except:
            print('\t- - - - - - - - - -')
            print(f'\tNo PanSTARRS data for {obj}.')
            return pd.DataFrame()
        
        else:
            lc.rename(columns = {'filterID' : 'filter', 'obsTime' : 'mjd', 'psfFlux' : 'flux_psf', 'psfFluxErr' : 'flux_psf_err', 'apFlux' : 'flux_aper', 'apFluxErr' : 'flux_aper_err'}, inplace = True)
            lc['filter'].replace({1 : 'g_ps', 2 : 'r_ps', 3 : 'i_ps', 4 : 'z_ps', 5 : 'y_ps'}, inplace = True)
            # data.query('abs(psfQfPerfect) > 0.5', inplace = True)
        
            for r, row in lc.iterrows():
                filt, flux_psf, flux_psf_err, flux_aper, flux_aper_err = row[['filter', 'flux_psf', 'flux_psf_err', 'flux_aper', 'flux_aper_err']]
                mag_psf = -2.5 * np.log10(flux_psf) + 8.90
                mag_psf_err = 1.08 * (flux_psf_err / flux_psf)
                mag_aper = -2.5 * np.log10(flux_aper) + 8.90
                mag_aper_err = 1.08 * (flux_aper_err / flux_aper)
                lc.loc[r, 'mag'] = mag_psf
                lc.loc[r, 'mag_err'] = mag_psf_err
                lc.loc[r, 'mag_aper'] = mag_aper
                lc.loc[r, 'mag_aper_err'] = mag_aper_err
            
            # for col in lc.columns:
            #     if col not in ['mjd', 'mag', 'mag_err', 'mag_aper', 'mag_aper_err', 'flux_psf', 'flux_psf_err', 'flux_aper', 'filter']:
            #         lc.drop(columns = col, inplace = True)
            
            lc['infoFlag'] -= (1+2+4+16+32+64+2097152+8388608+16777216+33554432+67108864+134217728)
            
            lc['flag'] = 0
            lc['source'] = 'panstarrs'
            lc.dropna(subset = ['mjd', 'mag'], inplace = True)
            lc = mag_to_flux(lc)
            lc.sort_values(by = 'mjd', inplace = True)
            lc.reset_index(drop = True, inplace = True)
                    
            N_filt = len(set(lc['filter']))
            print('\t- - - - - - - - - -')
            print(f'\tFound {len(lc.index)} Pan-STARRS photometry points in {N_filt} bands:')
            for filt in set(lc['filter']):
                lc_filt = lc.query('filter == @filt')
                N = len(lc_filt.index)
                fraction_flagged = np.sum(lc_filt['flag']) / N
                print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
            
            if save:
                lc.to_csv(savename, index = False)
            
            return lc