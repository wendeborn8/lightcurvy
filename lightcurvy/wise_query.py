import os
import pandas as pd
import numpy as np

from astroquery.ipac.irsa import Irsa
from astropy import units as u

from . import util
from .mag_to_flux import mag_to_flux

###

def wise_query(obj = None, ra = None, dec = None, search_radius = 0.0015,
                save = True, overwrite = False, obj_save = None, save_dir = None, **kwargs):
    
    print('Grabbing WISE light curve')
    
    if obj is not None and (ra is None and dec is None):
        ra, dec = util.get_radec(obj = obj)
        
    # Check if this light curve has already been generated
    if obj_save is None:
        obj_save = obj.lower().replace(' ', '_')
    if save_dir is None:
        save_dir = os.getcwd()
    savename = os.path.join(save_dir, 'wise', rf'{obj_save}.csv')
    os.makedirs(os.path.join(save_dir, 'wise'), exist_ok = True)
    
    if os.path.exists(savename) and not overwrite:
        lc = pd.read_csv(savename)
        print(f'WISE light curve already exists for {obj}, but will not overwrite. Grabbing light curve from {savename}.')
        return lc
    
    else:
        print(f'\tObject: {obj}\n\tEpochs: WISE, AllWISE, NeoWISE')
        lc_wise = Irsa.query_region(f'{ra} {dec}', catalog = 'allsky_4band_p1bs_psd', spatial = 'Cone', radius = search_radius*u.deg).to_pandas()
        lc_wise['source'] = 'wise'
        lc_wise.rename(columns = {'scan_id' : 'image_id'}, inplace = True)
        lc_allwise = Irsa.query_region(obj, catalog = 'allwise_p3as_mep', spatial = 'Cone', radius = search_radius*u.deg).to_pandas()
        lc_allwise['source'] = 'allwise'
        lc_allwise.rename(columns = {'w1mpro_ep' : 'w1mpro', 'w1sigmpro_ep' : 'w1sigmpro',
                                     'w1flux_ep' : 'w1flux', 'w1sigflux_ep' : 'w1sigflux',
                                     'w2mpro_ep' : 'w2mpro', 'w2sigmpro_ep' : 'w2sigmpro',
                                     'w2flux_ep' : 'w2flux', 'w2sigflux_ep' : 'w2sigflux',
                                     'frame_id' : 'image_id'}, inplace = True)
        lc_neowise = Irsa.query_region(obj, catalog = 'neowiser_p1bs_psd', spatial = 'Cone', radius = search_radius*u.deg).to_pandas()
        lc_neowise['source'] = 'neowise'
        lc_neowise.rename(columns = {'scan_id' : 'image_id'}, inplace = True)
    
        lc = pd.concat([lc_wise, lc_allwise, lc_neowise])
        
        for r, row in lc.iterrows():
            if row['cc_flags'] != '0000' or row['qual_frame'] == 0.0:
                flag = True
            else:
                flag = False
            lc.loc[r, 'flag'] = flag
        
        lcs_filts = []
        for col in ['w1mpro', 'w2mpro', 'w3mpro', 'w4mpro']:
            col_err = col.replace('mpro', 'sigmpro')
            filt = col[:2]
            lc_filt = lc.loc[:, ['mjd', 'source', 'image_id', 'flag', col, col_err]]
            lc_filt['filter'] = filt
            lc_filt.rename(columns = {col : 'mag', col_err : 'mag_err'}, inplace = True)
            lcs_filts.append(lc_filt)
        lc = pd.concat(lcs_filts)
        lc.dropna(subset = ['mag'], inplace = True)
        lc = mag_to_flux(lc)
        lc.sort_values(by = ['mjd'], inplace = True)
        lc.reset_index(drop = True, inplace = True)
        
        N_filt = len(set(lc['filter']))
        print('\t- - - - - - - - - -')
        print(f'\tFound {len(lc.index)} ZTF photometry points in {N_filt} bands:')
        for filt in set(lc['filter']):
            lc_filt = lc.query('filter == @filt')
            N = len(lc_filt.index)
            fraction_flagged = np.sum(lc_filt['flag']) / N
            print(f'\t\t{filt}: ', N, f'({fraction_flagged:.2f}% flagged)')
        print('\t\t- - - - - - - - - -')
        for source in set(lc['source']):
            lc_source = lc.query('source == @source')
            N = len(lc_source.index)
            fraction_flagged = np.sum(lc_source['flag']) / N
            print(f'\t\t{source}: ', N, f'({fraction_flagged:.2f}% flagged)')
        
        if save:
            lc.to_csv(savename, index = False)
    
        return lc