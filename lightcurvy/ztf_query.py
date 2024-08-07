import pandas as pd
import numpy as np
import os
import requests

from astropy.io.votable import parse

from . import util
from .mag_to_flux import mag_to_flux

###

def ztf_query(obj = None, ra = None, dec = None,
              save = True, overwrite = False, obj_save = None, save_dir = None,
              search_radius = 0.001,
              **kwargs):
    
    print('Grabbing ZTF light curve')
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'ztf', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'ztf'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'ZTF light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        url = fr'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE {ra} {dec} {search_radius}'
        print(f'\tObject: {obj}\n\tRA: {ra}\n\tDec: {dec}\n\tSearch Radius: {search_radius}\n\tZTF URL: {url}')
        
        response = requests.get(url)
        text = response.text
        with open('temp.txt', 'w+') as f:
            f.write(text)
        lc = parse('temp.txt').get_first_table().to_table(use_names_over_ids=True).to_pandas()
        
        if len(lc.index) == 0:
            print('\t- - - - - - - - - -')
            print(f'\tFound no ZTF photometry for {obj}.')
            return lc
        
        lc.rename(
            columns = {
                'oid' : 'object_id',
                # 'catflags' : 'flag',
                'filtercode' : 'filter',
                'magerr' : 'mag_err',
                'limitmag' : 'mag_limit',
                'expid' : 'image_id'}, 
            inplace = True
            )
        # lc.drop(columns = ['hjd', 'ra', 'dec', 'chi', 'sharp', 'filefracday', 'field', 'ccdid', 'qid', 'magzp', 'magzprms', 'clrcoeff', 'clrcounc', 'exptime', 'airmass', 'programid'], inplace = True)
        lc['flag'] = [0 if flag < 32768 else 1 for flag in lc['catflags']]
        lc['source'] = 'ztf'
        lc['filter'].replace({'zg' : 'g_ztf', 'zr' : 'r_ztf', 'i' : 'i_ztf'}, inplace = True)
        
        #
        lc = mag_to_flux(lc)        
        
        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} ZTF photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
        
        if save:
            lc.to_csv(savename, index = False)
        
        return lc