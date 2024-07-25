import io
import re
import sys
import time
import copy
import requests
import json
import os
import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from bs4 import BeautifulSoup as bs
import configparser

from astroquery.gaia import Gaia

###

# from lightcurvy import util
from lightcurvy.tess_query import tess_query
from lightcurvy.asassn_query import asassn_query
from lightcurvy.ztf_query import ztf_query
from lightcurvy.atlas_query import atlas_query
from lightcurvy.aavso_query import aavso_query
from lightcurvy.wise_query import wise_query
from lightcurvy.lcogt_query import lcogt_query
from lightcurvy.gaia_query import gaia_query
from lightcurvy.panstarrs_query import panstarrs_query
from lightcurvy.auto_period import auto_period

#%%

class Lightcurvy:
    
    def __init__(self, obj = None, obj_save = None, ra = None, dec = None, save_dir = None):
        
        self.obj = obj
        self.obj_save = obj_save
        self.ra = ra
        self.dec = dec
        self.save_dir = save_dir
        
        self._read_config()
            
    def query_all(self, 
                  save = True,
                  obj_save = None,
                  save_dir = None,
                  overwrite = False,
                  sources = ['asassn', 'ztf', 'atlas', 'aavso', 'wise', 'panstarrs', 'gaia', 'lcogt'], 
                  search_radius = 0.0015,
                  lcogt_start = None, lcogt_end = None,
                  overwrite_images = False,
                  tess_sectors = None, ffi_cutout_size = 50, limit_mag = 15, timeout = 6000,
                  ):
        
        self.save = save
        self.overwrite = overwrite
        self.search_radius = search_radius
        
        # 
        if obj_save is None and self.obj_save is None:
            self._make_obj_save()
            
        # Determine the paths
        self._make_dirs()
        
        
        # String to print between queries
        break_string = '#' * 50
        
        # Query each light curve source
        lcs, sources_successful = [], []
        for source, query_func in zip(['aavso', 'atlas', 'asassn', 'ztf', 'gaia', 'panstarrs', 'wise', 'lcogt'],
                                      [aavso_query, atlas_query, asassn_query, ztf_query, gaia_query, panstarrs_query, wise_query, lcogt_query]):
        # for source, query_func in zip(['asassn', 'ztf'],
        #                               [asassn_query, ztf_query]):
            
            if source not in sources:
                continue
            
            arguments = {
                'obj' : self.obj,
                'obj_save' : self.obj_save,
                'ra' : self.ra,
                'dec' : self.dec,
                'save' : self.save,
                'overwrite' : self.overwrite,
                'save_dir' : self.save_dir,
                'atlas_username' : self.atlas_username,
                'atlas_password' : self.atlas_password,
                'atlas_token' : self.atlas_token,
                'lcogt_username' : self.lcogt_username,
                'lcogt_password' : self.lcogt_password,
                'lcogt_token' : self.lcogt_token,
                'search_radius' : self.search_radius,
                'lcogt_start' : lcogt_start,
                'lcogt_end' : lcogt_end,
                'overwrite_images' : overwrite_images,
                'tess_sectors' : tess_sectors,
                'ffi_cutout_size' : ffi_cutout_size,
                'limit_mag' : limit_mag,
                'timeout' : timeout
                }
            
            try:
                lc = query_func(**arguments)
                
            except Exception as e:
                print(e)
                print(f'Issue retrieving {source} light curve for {self.obj}.')
                continue
                
            if len(lc.index) > 0:
                lcs.append(lc)
                sources_successful.append(source)
            
            print(break_string)
            
        # 
        print(f'Found light curve data for {len(sources_successful)} sources for {self.obj}:\n\t{sources_successful}')
        
        # Combine found light curves
        self.lc = pd.concat(lcs)
        
        # Restrict only to relevant columns
        self.lc = self.lc.loc[:, ['mjd', 'filter', 'source', 'flag', 'mag', 'mag_err', 'flux_Jy', 'err_Jy', 'flux_Flam', 'err_Flam', 'flux_lamFlam', 'err_lamFlam']]
        
        
        # Do some basic cleaning
        bad_filters = ['CV', 'H', 'HA', 'TG', 'CR', 'TB', 'TR', 'Vis.']
        self.lc.query('filter not in @bad_filters', inplace = True) # Drop irrelevant filters
        # lc.dropna(axis = 0, how = 'all', inplace = True) # Drop 
        self.lc.sort_values(by = ['mjd'], inplace = True)
        self.lc.reset_index(drop = True, inplace = True)
        
        # Drop anything 10 sigma above/below the median
        self.lc.query('mag > 0 and mag < 25', inplace = True)
        for filt in set(self.lc['filter']):
            lc_filt = self.lc.query('filter == @filt')
            med, std = np.nanmedian(lc_filt['mag']), np.nanstd(lc_filt['mag'])
            bad_inds = lc_filt.query('mag < @med-10*@std or mag > @med+10*@std').index
            self.lc.loc[bad_inds, 'flag'] = 1
            
        
        
    ###
    
    def multi_period(self, x_col = 'mjd', y_col = 'flux_Jy', ye_col = 'err_Jy', filters = None, 
                     period_range = [1.1, 14.9], mjd_min = None, mjd_max= None):
            
        lc = self.lc.copy()
        if mjd_min is not None:
            lc.query('mjd > @mjd_min', inplace = True)
        if mjd_max is not None:
            lc.query('mjd < @mjd_max', inplace = True)
            
        lc.query('flag != 1', inplace = True)
            
        if filters is None:
            filters = set(lc['filter'])
            
        fig, axs = plt.subplots(len(filters), 1, figsize = (6, len(filters) / 1.5), sharex = True, sharey = False)
            
        for filt, ax in zip(filters, axs):
            print(filt)
            lc_filt = self.lc.query('filter == @filt')
            if len(lc_filt.index) > 10:            
                results = auto_period(data = lc_filt, x_col = x_col, y_col = y_col, ye_col = ye_col, n_points = 2000, methods = ['ls'], period_range = period_range, n_mc = 100, n_passes = 1, plot = False, power_min = 0)
                
                period, power = results['Pass 1']['combined']['period'], results['Pass 1']['combined']['power']
                if len(results['Pass 1']['combined']['period_peaks']) > 0:
                    P = results['Pass 1']['combined']['period_peaks'][0]
                else:
                    P = np.nan
            else:
                period, power = np.nan, np.nan
                P = np.nan
                
            ax.plot(period, power, linewidth = 1)
            ax.axvline(P, linestyle = ':', color = 'red')
            ax.text(s = f'{len(lc_filt.index)} - {P:.2f}', x = 0.01, y = 0.95, ha = 'left', va = 'top', transform = ax.transAxes)
            ax.set_ylabel(filt)
            ax.set_xlim(min(period_range), max(period_range))
            ax.grid(alpha = 0.25)
        
        fig.subplots_adjust(hspace = 0.0)
        plt.show()
        
    ###
    
    def plot(self, x_col = 'mjd', y_col = 'flux_Jy', ye_col = 'err_Jy', filters = None, period_range = [0.1, 14.9], mjd_min = None, mjd_max= None, norm = False):
        
        filter_info = pd.read_excel(r'D:\My Drive\Data\filter_info.xlsx', index_col = 0)
        
        lc = self.lc.copy()
        if mjd_min is not None:
            lc.query('mjd > @mjd_min', inplace = True)
        if mjd_max is not None:
            lc.query('mjd < @mjd_max', inplace = True)
            
        lc.query('flag != 1', inplace = True)
            
        if filters is None:
            filters = set(lc['filter'])
        
        fig, ax = plt.subplots()
        for filt in filters:
            lc_filt = lc.query('filter == @filt')
            
            x = lc_filt[x_col]
            y = lc_filt[y_col]
            ye = lc_filt[ye_col]
            
            # 
            if norm:
                y /= np.nanmedian(y)
                ye /= np.nanmedian(y)
            
            ax.scatter(x, y, label = filt, s = 10, alpha = 0.75)
            
        ax.set_xlabel('Time')
        ax.set_ylabel('Value')
        
        ax.grid(alpha = 0.25)
        ax.legend(bbox_to_anchor = (1.0, 0.5), loc = 'center left')
        plt.show()
    
    ###
    
    def group_filters(self):
        
        filter_infos = pd.read_excel(r'D:\My Drive\Data\filter_info.xlsx', index_col = 0)
        groups = ['u', 'B', 'g', 'V', 'r', 'i', 'z', 'y_ps', 'w1']
        
        # for group in group:
        
        self.lc['filter_orig'] = list(self.lc['filter'])
        self.lc['filter'] = np.array([filter_infos.loc[filt, 'filter_group'] for filt in self.lc['filter_orig']])
            
    ###
    
    def _mag_to_flux(self):
        from mag_to_flux import mag_to_flux
        self.lc = mag_to_flux(self.lc)
        
    def _make_obj_save(self):
        bad_chars = ['*', '/', '\\', ':', '?', '\"', '<', '>', '|', '.', ' ']
        replace_with = '_'
        
        self.obj_save = copy.deepcopy(self.obj)
        for bad_char in bad_chars:
            self.obj_save.replace(bad_char, replace_with)
        
        return
    
    def _make_dirs(self, save_dir = None):
        # save_dir
            # save_dir/lightcurves
            # save_dir/*source*
        if self.save_dir is None:
            if save_dir is None:
                self.save_dir = r'./' # Current directory
            else:
                self.save_dir = save_dir
            
        self.lightcurves_dir = os.path.join(self.save_dir, 'lightcurves')
        
    def _read_config(self):
        config_path = r'D:\My Drive\Python\lightcurvy\config.txt'
        config = configparser.ConfigParser()
        config.read(config_path)
        
        self.atlas_username = config.get('Atlas', 'atlas_username')
        self.atlas_password = config.get('Atlas', 'atlas_password')
        self.atlas_token = config.get('Atlas', 'atlas_token')
        self.lcogt_username = config.get('LCOGT', 'lcogt_username')
        self.lcogt_password = config.get('LCOGT', 'lcogt_password')
        self.lcogt_token = config.get('LCOGT', 'lcogt_token')
        
        