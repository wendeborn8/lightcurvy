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
from .tess_query import tess_query
from .asassn_query import asassn_query
from .ztf_query import ztf_query
from .atlas_query import atlas_query
from .aavso_query import aavso_query
from .wise_query import wise_query
from .lcogt_query import lcogt_query
from .gaia_query import gaia_query
from .panstarrs_query import panstarrs_query
from .auto_period import auto_period

###

class Lightcurvy:
    
    def __init__(self, obj = None, obj_save = None, ra = None, dec = None, save = None, overwrite = None, save_dir = None):
        
        '''
        
        
        Paramaters
        ----------
        obj : str
            Object name to use when searching for light curves. Cross-matches SIMBAD, so make sure to provide SIMBAD compliant name. Currently, you MUST proivde 'obj'.
            
        obj_save : str
            Prefix for all file names to use when saving files. If not proivded, will replace special characters and spaces in 'obj' with '_'. 
            
        ra : float
            Right Ascension of your target.
            
        dec : float
            Declination of your target.
            
        save : bool, default True
            Whether to save the resulting individual and combined light curves.
            
        save_dir : str
            Where to save the light curves. Sub-directories will be created for individual sources. If not provided, will use the current working directory
            
        overwrite : bool, default False
            Whether to overwrite any light curves that already exist. Does not control overwriting LCOGT images or intermediate files for TESS ligth curves, will only redo creation of light curves.
            
        
        Returns
        ----------
        None
            
        '''
        
        self.obj = obj
        self.obj_save = obj_save
        self.ra = ra
        self.dec = dec
        self.save = save
        self.overwrite = overwrite
        self.save_dir = save_dir
        
        self._read_config()
        self._read_filter_infos()
            
    def query_all(self, 
                  save = True,
                  obj_save = None,
                  save_dir = None,
                  overwrite = False,
                  sources = ['asassn', 'ztf', 'atlas', 'aavso', 'wise', 'panstarrs', 'gaia', 'lcogt', 'tess'], 
                  search_radius = 0.0015,
                  lcogt_start = None, lcogt_end = None,
                  overwrite_images = False,
                  tess_sources = ['tglc', 'spoc', 'eleanor'], 
                  tess_sectors = None, ffi_cutout_size = 50, limit_mag = 15, timeout = 6000,
                  eleanor_size = 15, eleanor_bkg_size = 31
                  ):
        
        '''
        Method to query all sources (or at least those specified) for a target
        
        Parameters
        ----------
        save : bool, default True
            Whether to save the resulting individual and combined light curves.
            
        obj_save : str
            Prefix for all file names to use when saving files. If not proivded, will replace special characters and spaces in 'obj' with '_'. 
            
        save_dir : str
            Where to save the light curves. Sub-directories will be created for individual sources. If not provided, will use the current working directory
            
        overwrite : bool, default False
            Whether to overwrite any light curves that already exist. Does not control overwriting LCOGT images or intermediate files for TESS ligth curves, will only redo creation of light curves
            
        sources : list of str, default ['asassn', 'ztf', 'atlas', 'aavso', 'wise', 'panstarrs', 'gaia', 'lcogt', 'tess']
            What light curve sources to query. Options include: 'aavso', 'asassn', 'atlas', 'ztf', 'wise', 'panstarrs', 'tess', 'lcogt', and 'gaia'.
            
        search_radius : float, default 0.0015
            Radius (in degrees) to use for certain queries when searching for data. Default is 0.0015 degrees (5.4 arcseconds)
            
        lcogt_start : str
            Start date for LCOGT query, in format 'YYYY-MM-DD'. If not provided, default to '2014-01-01'
            
        lcogt_end : str
            Start date for LCOGT query, in format 'YYYY-MM-DD'. If not provided, default to the most recent date.
            
        overwrite_images : bool, default False
            Whether to overwrite images download from LCOGT. Recommended to set to True ONLY when a substantuce change has been made to LCOGT's BANZAI pipeline that affects archived frames.
            
        tess_sources : list of str, default ['tglc', 'spoc', 'eleanor']
            What TESS light curve sources to query. Options include 'spoc', 'eleanor', 'tglc'. If you only need a quick light curve for older sectors, recommend using only 'spoc'. Technical issues often arise with 'eleanor' and 'tglc' can be very slow and not always reliable.
            
        tess_sectors : list of int
            What TESS sectors to fetch light curves for. If not provided, will retrieve all possible sectors.
        
        ffi_cutout_size : int, default 50
            What size cutout to use when running 'tglc'. 50+ is recommended for best results.
            
        limit_mag : float, default 16
            What Gaia G magnitude to limit background source subtraction to when running 'tglc'. Larger values (i.e. dimmer stars) should result in better results, but will take longer.
            
        timeout : int, default 6000
            ### NOT CURRENTLY IMPLEMENTED! Must manually change timeout parameter in XX file in Astroquery ### 
            Controls the MAST Astroquery timeout, in seconds. By default it is 600 seconds, which isn't enough for 'tglc' to fetch the FFIs. 
            
        eleanor_size : int, default 15
            Size to use for running 'eleanor'. Must be odd integer. Same as 'width' and 'height'. See 'eleanor' documentation for more info.
            
        eleanor_bkg_size : int, default 31
            Bakcground size to use for running 'eleanor'. Must be odd integer. Same as 'bkg_size'. See 'eleanor' documentation for more info.
            
            
        Returns
        ----------
        DataFrame
            Combined light curve as a Pandas DataFrame
        
        '''
        
        if self.save is None:
            self.save = save
        if self.overwrite is None:
            self.overwrite = overwrite
        self.search_radius = search_radius
        
        # 
        if obj_save is None and self.obj_save is None:
            self._make_obj_save()
            
        # Determine the paths
        self._make_dirs()
        
        
        # String to print between queries
        break_string = '#' * 50
        
        source_query_funcs = {'aavso' : aavso_query, 'asassn' : asassn_query, 'atlas' : atlas_query, 'ztf' : ztf_query, 'gaia' : gaia_query, 'panstarrs' : panstarrs_query, 'wise' : wise_query, 'lcogt' : lcogt_query, 'tess' : tess_query}
        # source_query_funcs = {'ztf' : ztf_query}
        
        # Query each light curve source
        lcs, sources_successful = [], []
        for source in sources:
            if source not in source_query_funcs.keys():
                print(f'Source {source} not recognized. Acceptable sources include: {source_query_funcs.keys}.')
                continue
            
            query_func = source_query_funcs[source]
            
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
                print(break_string)
                continue
                
            if len(lc.index) > 0:
                lcs.append(lc)
                sources_successful.append(source)
            
            print(break_string)
            
        # 
        print(f'Found light curve data for {len(sources_successful)} sources for {self.obj}:\n\t{sources_successful}')
        
        # Combine found light curves
        self.lc = pd.concat(lcs)
        
        # Do some basic cleaning
        bad_filters = ['CV', 'H', 'HA', 'TG', 'CR', 'TB', 'TR', 'Vis.']
        self.lc.query('filter not in @bad_filters', inplace = True) # Drop irrelevant filters
        self.lc.sort_values(by = ['mjd'], inplace = True)
        self.lc.reset_index(drop = True, inplace = True)
        
        # Drop anything 10 sigma above/below the median
        self.lc.query('mag > 0 and mag < 25', inplace = True)
        for filt in set(self.lc['filter']):
            lc_filt = self.lc.query('filter == @filt')
            med, std = np.nanmedian(lc_filt['mag']), np.nanstd(lc_filt['mag'])
            bad_inds = lc_filt.query('mag < @med-10*@std or mag > @med+10*@std').index
            self.lc.loc[bad_inds, 'flag'] = 1
            
        #
        if self.save:
            self.lc.to_csv(os.path.join(self.save_dir, 'lightcurves', f'{self.obj_save}.csv'), index = False)
            
        return self.lc
            
        
        
    ###
    
    def multi_period(self, x_col = 'mjd', y_col = 'flux_Jy', ye_col = 'err_Jy', filters = None, 
                     period_range = [0.1, 14.9], mjd_min = None, mjd_max= None, flag = False):
        
        '''
        Perform a rough periodogram for all the filters in the resulting light curve
        
        Parameters
        ----------
        x_col : str, default 'mjd'
            Column name to use for the time axis.
            
        y_col : str, default 'flux_Jy'
            Column name to use for the y values in the periodogram.
            
        ye_col : str, default 'err_Jy'
            Column name to use for the y value uncertainties in the periodogram. Only utilized if 'n_mc' is greater than 1
            
        filters : list of str
            Which filters to use for the periodogram. If not provided, will use all filters.
            
        period_range : list/array of length 2, default [0.1, 14.9]
            What range of periods to search for signals. Anything larger than ~1/2 the total time baseline is typically not useful, unless the signal is VERY coherent.
            
        mjd_min : float
            Start MJD to cut off the data at.
            
        mjd_max : float
            End MJD to cut off the data at.
            
        flag : bool, default False
            Whether to remove flagged data. True: Remove flagged data, false: use all data.
            
            
        Returns
        ----------
        None
            
        '''
        
        lc = copy.deepcopy(self.lc)
        if mjd_min is not None:
            lc.query('mjd > @mjd_min', inplace = True)
        if mjd_max is not None:
            lc.query('mjd < @mjd_max', inplace = True)
        
        if flag:
            lc.query('flag != 1', inplace = True)
            
        if filters is None:
            filters = set(lc['filter'])
            
        fig, axs = plt.subplots(len(filters), 1, figsize = (6, len(filters)), sharex = True, sharey = False)
            
        for filt, ax in zip(filters, axs):
            print(filt)
            filt_info = self.filter_infos.loc[filt, :]
            color = filt_info['color']
            lc_filt = lc.query('filter == @filt')
            if len(lc_filt.index) > 10:            
                results = auto_period(data = lc_filt, x_col = x_col, y_col = y_col, ye_col = ye_col, n_points = 1000, methods = ['ls', 'pdm'], period_range = period_range, n_mc = 100, n_passes = 1, rand_frac = 0.75, plot = False, power_min = 0)
                
                period, power = results['Pass 1']['combined']['period'], results['Pass 1']['combined']['power']
                if len(results['Pass 1']['combined']['period_peaks']) > 0:
                    P = results['Pass 1']['combined']['period_peaks'][0]
                else:
                    P = np.nan
            else:
                period, power = np.nan, np.nan
                P = np.nan
                
            ax.plot(period, power, linewidth = 1, color = color)
            ax.axvline(P, linestyle = ':', color = 'red')
            ax.text(s = f'{len(lc_filt.index)} - {P:.2f}', x = 0.01, y = 0.95, ha = 'left', va = 'top', transform = ax.transAxes)
            ax.set_ylabel(filt)
            ax.set_xlim(min(period_range), max(period_range))
            ax.grid(alpha = 0.25)
        
        fig.subplots_adjust(hspace = 0.0)
        plt.show()
        
    ###
    
    def plot_lc(self, x_col = 'mjd', y_col = 'flux_Jy', ye_col = 'err_Jy', filters = None, mjd_min = None, mjd_max= None, norm = False, flag = False):
        
        '''
        Plot the resulting light curve, separating each filter
        
        Paramaters
        ----------
        x_unit : str, default 'um'
            What unit to use for the x axis. Options include 'um', 'A', 'nm', 'GHz'.
        
        y_col : str, default 'flux_Jy'
            Column name to use for the y values.
            
        ye_col : str, default 'err_Jy'
            Column name to use for the y value uncertainties.
        
        filters : list of str
            Which filters to use for the periodogram. If not provided, will use all filters.
        
        mjd_min : float
            Start MJD to cut off the data at.
            
        mjd_max : float
            End MJD to cut off the data at.
            
        flag : bool, default False
            Whether to remove flagged data. True: Remove flagged data, false: use all data.
            
        norm : bool, default False
            Whether to normalize each filter by the median. Can sometimes make for easier comparison.
            
            
        Returns
        ----------
        None
        
        '''
        
        lc = self.lc.copy()
        if mjd_min is not None:
            lc.query('mjd > @mjd_min', inplace = True)
        if mjd_max is not None:
            lc.query('mjd < @mjd_max', inplace = True)
        
        if flag:
            lc.query('flag != 1', inplace = True)
            
        if filters is None:
            filters = set(lc['filter'])
        
        fig, ax = plt.subplots()
        for filt in filters:
            filt_info = self.filter_infos.loc[filt, :]
            color, marker = filt_info['color'], filt_info['marker']
            
            lc_filt = lc.query('filter == @filt')
            
            x = lc_filt[x_col]
            y = lc_filt[y_col]
            ye = lc_filt[ye_col]
            
            # 
            if norm:
                y /= np.nanmedian(y)
                ye /= np.nanmedian(y)
            
            ax.scatter(x, y, label = filt, s = 10, alpha = 0.75, edgecolor = color, facecolor = 'None', marker = marker)
            
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        if y_col == 'mag':
            ax.invert_yaxis()
            
        ax.text(s = self.obj, x = 0.02, y = 0.95, ha = 'left', va = 'top', transform = ax.transAxes)
        ax.tick_params(direction = 'in', right = True, top = True)
        ax.grid(alpha = 0.25)
        ax.legend(bbox_to_anchor = (0.5, 1.0), loc = 'lower center', ncol = 6)
        plt.show()
        
    ###
    
    def plot_sed(self, x_unit = 'um', y_col = 'flux_Flam', ye_col = 'err_Flam', filters = None, mjd_min = None, mjd_max= None, flag = False):
        
        '''
        Plot an SED of the resulting light curve, separating each filter
        
        Paramaters
        ----------
        x_unit : str, default 'um'
            What unit to use for the x axis. Options include 'um', 'A', 'nm', 'GHz'.
        
        y_col : str, default 'flux_Jy'
            Column name to use for the y values.
            
        ye_col : str, default 'err_Jy'
            Column name to use for the y value uncertainties.
        
        filters : list of str
            Which filters to use for the periodogram. If not provided, will use all filters.
        
        mjd_min : float
            Start MJD to cut off the data at.
            
        mjd_max : float
            End MJD to cut off the data at.
            
        flag : bool, default False
            Whether to remove flagged data. True: Remove flagged data, false: use all data.
            
            
        Returns
        ----------
        None
        
        '''
        
        
        lc = self.lc.copy()
        if mjd_min is not None:
            lc.query('mjd > @mjd_min', inplace = True)
        if mjd_max is not None:
            lc.query('mjd < @mjd_max', inplace = True)
        
        if flag:
            lc.query('flag != 1', inplace = True)
            
        if filters is None:
            filters = set(lc['filter'])
        
        fig, ax = plt.subplots()
        for filt in filters:
            filt_info = self.filter_infos.loc[filt, :]
            wl_eff, color, marker = filt_info['wl_eff'], filt_info['color'], filt_info['marker']
            
            lc_filt = lc.query('filter == @filt')
            
            y = lc_filt[y_col]
            ye = lc_filt[ye_col]
            
            if x_unit == 'um':
                wl_eff /= 1e4
            if x_unit == 'nm':
                wl_eff /= 1e3
            if x_unit == 'GHz':
                wl_eff  = (1 / wl_eff) * (1e10) * (3e8) * (1e-9)
                
            x = np.array([wl_eff for _ in y])
            
            ax.scatter(x, y, s = 5, alpha = 0.2, edgecolor = 'None', facecolor = color, marker = marker)
            ax.errorbar(wl_eff, np.nanmedian(y), yerr = np.nanstd(y), linestyle = '', ecolor = color, capsize = 2,
                       label = filt, markersize = 7.5, alpha = 1, markeredgecolor = color, markerfacecolor = 'white', marker = marker)
        
        if x_unit in ['A', 'nm', 'um']:
            x_label = f'Wavelength ({x_unit})'
        else:
            x_label = f'Frequency ({x_unit})'
        
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_col)
        if y_col == 'mag':
            ax.invert_yaxis()
        
        ax.tick_params(direction = 'in', right = True, top = True, which = 'both')
        ax.text(s = self.obj, x = 0.98, y = 0.95, ha = 'right', va = 'top', transform = ax.transAxes)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(alpha = 0.25)
        ax.legend(bbox_to_anchor = (0.5, 1.0), loc = 'lower center', ncol = 6)
        plt.show()
    
    ###
    
    def group_filters(self):
        
        groups = ['u', 'B', 'g', 'V', 'r', 'i', 'z', 'y_ps', 'w1']
                
        self.lc['filter_orig'] = list(self.lc['filter'])
        self.lc['filter'] = np.array([self.filter_infos.loc[filt, 'filter_group'] for filt in self.lc['filter_orig']])
            
    ###
        
    def _make_obj_save(self):
        bad_chars = ['*', '/', '\\', ':', '?', '\"', '<', '>', '|', '.', ' ']
        replace_with = '_'
        
        self.obj_save = copy.deepcopy(self.obj)
        for bad_char in bad_chars:
            self.obj_save.replace(bad_char, replace_with)
            
    ###
            
    def _make_dirs(self, save_dir = None):
        if self.save_dir is None:
            if save_dir is None:
                self.save_dir = r'./' # Current directory
            else:
                self.save_dir = save_dir
            
        self.lightcurves_dir = os.path.join(self.save_dir, 'lightcurves')
        
    ###
        
    def modify_config(self, atlas_username = None, atlas_password = None, atlas_token = None, lcogt_username = None, lcogt_password = None, lcogt_token = None):
        config_path = os.path.join(os.path.dirname(__file__), 'config.txt')
        config = configparser.ConfigParser()
        config.read(config_path)
        
        if atlas_username: config['Atlas']['atlas_username'] = atlas_username
        if atlas_password: config['Atlas']['atlas_password'] = atlas_password
        if atlas_token: config['Atlas']['atlas_token'] = atlas_token
        if lcogt_username: config['LCOGT']['lcogt_username'] = lcogt_username
        if lcogt_password: config['LCOGT']['lcogt_password'] = lcogt_password
        if lcogt_token: config['LCOGT']['lcogt_token'] = lcogt_token
        
        with open(config_path, 'w') as configfile:
            config.write(configfile)
            
        self._read_config()
        
    ###
        
    def _read_config(self):
        config_path = os.path.join(os.path.dirname(__file__), 'config.txt')
        print(config_path)
        config = configparser.ConfigParser()
        config.read(config_path)
        
        self.atlas_username = config.get('Atlas', 'atlas_username')
        self.atlas_password = config.get('Atlas', 'atlas_password')
        self.atlas_token = config.get('Atlas', 'atlas_token')
        self.lcogt_username = config.get('LCOGT', 'lcogt_username')
        self.lcogt_password = config.get('LCOGT', 'lcogt_password')
        self.lcogt_token = config.get('LCOGT', 'lcogt_token')
        
    ###
        
    def _read_filter_infos(self):
        self.filter_infos_file = os.path.join(os.path.dirname(__file__), 'filter_infos.csv')
        self.filter_infos = pd.read_csv(self.filter_infos_file, index_col = 0)
        
        
