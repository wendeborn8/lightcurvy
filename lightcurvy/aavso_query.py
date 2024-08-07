import io
import re
import sys
import time
import requests
import os
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import tqdm
from bs4 import BeautifulSoup as bs

from .mag_to_flux import mag_to_flux

###

def aavso_query(obj = None, ra = None, dec = None, save_dir = None,
                save = True, overwrite = False, obj_save = None,
                **kwargs):
    
    print('Grabbing AAVSO light curve')
    
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'aavso', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'aavso'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'AAVSO light curve already exists for {obj} ({savename}), but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        if 'V*' in obj:
            obj_url = obj.replace('V*', '').strip().replace(' ', '+')
        else:
            obj_url = obj.replace(' ', '+')
        
        page_number = 1
        per_page = 200
        n_page_max = 100
        dfs = []
        mjd_min_prev = 99999
        mjd_min_abs = 50000
        
        print(f'\tObject: {obj}\n\tObs. Types: CCD\n\tMJD Min: {mjd_min_abs}')
        print('\t- - - - - - - - - -')
        
        pbar = tqdm(desc = '\tIterating through pages', total = n_page_max)
        while page_number < n_page_max:
            url = f'https://app.aavso.org/webobs/results/?star={obj_url}&num_results={per_page}&obs_types=ccd&page={page_number}'
            response = requests.get(url)
            stuff = response.text
            try:
                soup = bs(stuff, features = 'lxml')
            except:
                soup = bs(stuff, features = 'html.parser')
            rows_even = soup.find_all(attrs = {'class' : 'obs tr-even'})
            rows_odd = soup.find_all(attrs = {'class' : 'obs tr-odd'})
            rows = rows_odd + rows_even
            df = pd.DataFrame()
            for r, row in enumerate(rows):
                text = row.text
                entries = text.split('\n')
                
                aavso_obj, jd, mag, mag_err, filt, observer = entries[2], entries[3], entries[5], entries[6], entries[7], entries[8]
                try:
                    mag = float(mag)
                except:
                    mag = np.nan
                try:
                    mag_err = float(mag_err)
                except:
                    mag_err = np.nan
                index = f'{page_number}_{jd}'
                df.loc[index, ['obj', 'mjd', 'mag', 'mag_err', 'filter', 'observer']] = aavso_obj, float(jd)-2400000.5, float(mag), float(mag_err), filt, observer
            page_number += 1
            pbar.update(1)
            if len(df.index) == 0:
                break
            
            dfs.append(df)
            mjd_min = min(df['mjd'])

            if mjd_min < mjd_min_abs:
                break
            elif mjd_min == mjd_min_prev:
                break
            mjd_min_prev = copy.deepcopy(mjd_min)
            
        try: 
            lc = pd.concat(dfs)
        except:
            print('\t- - - - - - - - - -')
            print(f'\tFound no AAVSO photometry for {obj}')
            lc = pd.DataFrame()
            if save:
                lc.to_csv(savename, index = False)
            return lc
        
        # Replace some filter names, remove other filters entirely
        lc['filter'].replace({'SU' : 'up', 'SG' : 'gp', 'SR' : 'rp', 'SI' : 'ip', 'SZ' : 'zp',
                              'ZS' : 'zs', 
                              'STU' : 'u_stromgren', 'STV' : 'v_stromgren', 'STB' : 'b_stromgren', 'STY' : 'y_stromgren'}, inplace = True)
        bad_filters = ['CV', 'H', 'HA', 'TG', 'CR', 'TB', 'TR', 'Vis.']
        lc.query('filter not in @bad_filters', inplace = True)
        
        if len(lc.index) == 0:
            print('\n\t- - - - - - - - - -')
            print(f'\tNo AAVSO data found for {obj}')
            lc = pd.DataFrame()
            if save:
                lc.to_csv(savename, index = False)
            return lc
        
        lc['source'] = 'aavso'
        lc['flag'] = 0
        lc.dropna(subset = ['mjd', 'mag'], inplace = True)
        lc = mag_to_flux(lc)
        lc.sort_values(by = ['mjd'], inplace = True)
        lc.reset_index(inplace = True, drop = True)
                
        # Group very nearby observations
        lcs_binned = []
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt').copy()
            if len(lc_filt.index) > 3:
            
                mjd_min_filt, mjd_max_filt = min(lc_filt['mjd']), max(lc_filt['mjd'])
                bins = np.arange(mjd_min_filt, mjd_max_filt, 60/24/60) # 10 minute bins
                
                bin_med, bin_edges, binnumber = stats.binned_statistic(np.array(lc_filt['mjd']).astype(float), np.array(lc_filt['mag']).astype(float), statistic=np.nanmedian, bins = bins)
                bin_err, _, _ = stats.binned_statistic(np.array(lc_filt['mjd']).astype(float), np.array(lc_filt['mag']).astype(float), statistic=np.median, bins = bins)
                bin_centers = bin_edges[1:] - (bin_edges[1] - bin_edges[0])/2
                
                lc_binned = pd.DataFrame({'mjd' : bin_centers, 'filter' : filt, 'source' : 'aavso', 'flag' : 0, 'mag' : bin_med, 'mag_err' : bin_err, 'bin_width' : 1/24})
                lcs_binned.append(lc_binned)
                

        lc = pd.concat(lcs_binned)
        lc.dropna(subset = ['mjd', 'mag'], inplace = True)
        lc = mag_to_flux(lc)
        lc.sort_values(by = ['mjd'], inplace = True)
        lc.reset_index(inplace = True, drop = True)
        
        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} AAVSO photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.nan #np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
        
        if save:
            lc.to_csv(savename, index = False)
        
        return lc